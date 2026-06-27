#include "benchmark.h"
#include "estruturas.h"
#include "util.h"
#include "comum.h"
#include <cstdio>
#include <chrono>
#include <vector>
#include <random>
#include <algorithm>
#include <fstream>
using namespace std;

// =============================================================================
//  11. BENCHMARK COMPLETO
//  Metodologia: (1) so compara estruturas comparaveis (Grupo A = busca por
//  chave; Bloom no Grupo B); (2) alvo variado entre 64 pontos (evita O(1)
//  enganoso e hoisting); (3) repeticoes adaptativas p/ superar o timer (~0.6ms).
// =============================================================================

// Exibe uma linha de benchmark com barra proporcional -usa nanosegundos internamente
static void barraComparativa(const char *nome, double ns, double mx, int larg = 28)
{
    printf("  %-28s ", nome);
    printBarra(ns, mx, larg);
    if (ns < 1000.0)
        printf(" %8.2f ns/op\n", ns);
    else
        printf(" %8.2f us/op\n", ns / 1000.0);
}

// --- Medicao robusta de remocao em hash ---
static NoHash *g_br_baldes[TAMANHO_TABELA];
static HashSlot g_br_slots[TAM_HASH_OTM];

static double medirRemChain(const vector<unsigned int> &ids, const vector<double> &ar)
{
    int n = (int)ids.size();
    for (int b = 0; b < TAMANHO_TABELA; b++)
        g_br_baldes[b] = NULL;
    for (int i = 0; i < n; i++)
    {
        unsigned int id = ids[i];
        int b = id % TAMANHO_TABELA;
        NoHash *x = new NoHash;
        x->sample_id = id;
        x->bounding_area = ar[i];
        x->prox = g_br_baldes[b];
        g_br_baldes[b] = x;
    }
    asm volatile("" ::: "memory");
    auto t1 = chrono::high_resolution_clock::now();
    for (int i = 0; i < n; i++)
    {
        unsigned int id = ids[i];
        int b = id % TAMANHO_TABELA;
        NoHash *p = g_br_baldes[b], *ant = NULL;
        while (p && p->sample_id != id)
        {
            ant = p;
            p = p->prox;
        }
        if (p)
        {
            if (!ant)
                g_br_baldes[b] = p->prox;
            else
                ant->prox = p->prox;
            delete p;
        }
    }
    long long tot = chrono::duration_cast<chrono::nanoseconds>(
                        chrono::high_resolution_clock::now() - t1)
                        .count();
    return (double)tot / n;
}

static double medirRemOpen(const vector<unsigned int> &ids, const vector<double> &ar)
{
    int n = (int)ids.size();
    for (int s = 0; s < TAM_HASH_OTM; s++)
    {
        g_br_slots[s].ocupado = false;
        g_br_slots[s].removido = false;
    }
    for (int i = 0; i < n; i++)
    {
        unsigned int id = ids[i];
        int s = id % TAM_HASH_OTM;
        while (g_br_slots[s].ocupado)
            s = (s + 1) % TAM_HASH_OTM;
        g_br_slots[s] = {id, ar[i], true, false};
    }
    asm volatile("" ::: "memory");
    auto t1 = chrono::high_resolution_clock::now();
    for (int i = 0; i < n; i++)
    {
        unsigned int id = ids[i];
        int s = id % TAM_HASH_OTM, ini = s;
        while (g_br_slots[s].ocupado || g_br_slots[s].removido)
        {
            if (g_br_slots[s].ocupado && g_br_slots[s].sample_id == id)
            {
                g_br_slots[s].ocupado = false;
                g_br_slots[s].removido = true;
                break;
            }
            s = (s + 1) % TAM_HASH_OTM;
            if (s == ini)
                break;
        }
    }
    long long tot = chrono::duration_cast<chrono::nanoseconds>(
                        chrono::high_resolution_clock::now() - t1)
                        .count();
    // Varre o array inteiro contando tombstones -> observa TODAS as remocoes
    long long chk = 0;
    for (int s = 0; s < TAM_HASH_OTM; s++)
        if (g_br_slots[s].removido)
            chk++;
    g_sink_d += chk;
    return (double)tot / n;
}

void executarBenchmark(NoLista *lista, TabelaHashOtm *hotm, NoAVL *avl,
                       TabelaHash *hash, BloomFilter *bloom, unsigned int id_ref)
{
    // Conta a lista e usa o elemento do meio como ID de referencia no cabecalho.
    int n_lista = 0;
    for (NoLista *p = lista; p; p = p->prox)
        n_lista++;
    NoLista *meio = lista;
    for (int i = 0; i < n_lista / 2 && meio->prox; i++)
        meio = meio->prox;
    id_ref = meio->sample_id;

    // Coleta 64 alvos (id, area) distribuidos pelo dataset. Variar o alvo a cada
    const int NP = 64;
    unsigned int pid[NP];
    double parea[NP];
    {
        NoLista *p = lista;
        int passo = max(1, n_lista / NP);
        for (int k = 0; k < NP; k++)
        {
            // avanca 'passo' nos (ou reinicia se a lista acabar)
            for (int s = 0; s < passo && p && p->prox; s++)
                p = p->prox;
            if (!p)
                p = lista;
            pid[k] = p->sample_id;
            parea[k] = p->bounding_area;
        }
    }

    // Mede ns/op com reps ADAPTATIVO. O alvo cicla via volatile (sem hoisting).
    auto medirR = [&](int reps, auto fn) -> double
    {
        volatile int ctr = 0;
        auto t1 = chrono::high_resolution_clock::now();
        for (int i = 0; i < reps; i++)
        {
            int k = (ctr++) % NP;
            fn(pid[k], parea[k]);
        }
        long long tot = chrono::duration_cast<chrono::nanoseconds>(
                            chrono::high_resolution_clock::now() - t1)
                            .count();
        return (double)tot / reps;
    };
    const int REP_LENTO = 3000, REP_RAPIDO = 2000000;
    // cronoNs(): cronometra um bloco e divide pelo nº de operacoes -> ns/op.
    auto cronoNs = [&](int divisor, auto fn) -> double
    {
        auto t1 = chrono::high_resolution_clock::now();
        fn();
        asm volatile("" ::: "memory");
        return (double)chrono::duration_cast<chrono::nanoseconds>(
                   chrono::high_resolution_clock::now() - t1)
                   .count() /
               divisor;
    };

    double ns_henc = medirR(REP_RAPIDO, [&](unsigned int id, double)
                            { double a=0; g_sink_b = buscarHash(hash, id, &a); g_sink_d=a; });
    double ns_hotm = medirR(REP_RAPIDO, [&](unsigned int id, double)
                            { double a=0; g_sink_b = buscarHashOtm(hotm, id, &a); g_sink_d=a; });
    double ns_avl = medirR(REP_RAPIDO, [&](unsigned int, double ar)
                           { g_sink_b = (buscarAVL(avl, ar) != NULL); });
    double ns_lst = medirR(REP_LENTO, [&](unsigned int, double ar)
                           { g_sink_b = buscarLista(lista, ar); });
    double ns_blm = medirR(REP_RAPIDO, [&](unsigned int id, double)
                           { g_sink_b = buscarBloom(bloom, id); });

    printf("\n" CYN BOLD "  BENCHMARK COMPLETO - todas as estruturas implementadas\n" RST);
    printf(DIM "  (cada estrutura medida isoladamente; valores em ns por operacao)\n" RST);
    printLinha('=', 57);

    // GRUPO A: estruturas que respondem a MESMA pergunta ("achar chave X").
    double mx_a = max({ns_henc, ns_hotm, ns_avl, ns_lst, 0.01});
    printf(BOLD "\n  [ GRUPO A - BUSCA POR CHAVE (comparavel) ]" RST);
    printf(DIM " alvo: ID %u\n" RST, id_ref);
    printLinha('-', 57);
    barraComparativa("Hash Encadeamento O(1):", ns_henc, mx_a);
    barraComparativa("Hash Open Addr. (otim.):", ns_hotm, mx_a);
    barraComparativa("AVL O(log n):", ns_avl, mx_a);
    barraComparativa("Lista O(n) [R21 degradado]:", ns_lst, mx_a);
    printf(DIM "  Todas resolvem 'encontrar o elemento de chave X' -> comparacao valida.\n" RST);

    // GRUPO B: estrutura de PROPOSITO DISTINTO. NAO comparar na mesma regua.
    printf(BOLD "\n  [ GRUPO B - TESTE DE PERTINENCIA (NAO comparavel) ]\n" RST);
    printLinha('-', 57);
    barraComparativa("Bloom Filter O(k):", ns_blm, max(ns_blm, 0.01));
    printf(DIM "  Bloom responde 'X PODE estar no conjunto?' (com falsos positivos),\n"
               "  nao 'onde esta X?'. Por isso NAO entra no ranking do Grupo A.\n"
               "  Idem KD-Tree (busca espacial 2D) - ver opcao 2.\n" RST);

    // ===== INSERCAO: TODAS as estruturas (em copias TEMPORARIAS isoladas) =====
    const int N_INS = 300000;
    vector<double> areas(N_INS);
    {
        mt19937 g(99);
        uniform_real_distribution<double> d(0.001, 0.999);
        for (double &a : areas)
            a = d(g);
    }

    NoLista *t_lst = NULL;
    double ins_lst = cronoNs(N_INS, [&]
                             { for (int i=0;i<N_INS;i++) t_lst = inserirLista(t_lst, areas[i], (unsigned)i); });
    liberarLista(t_lst);

    double ins_avl = 0;
    {
        EstadoAVLGuard g;
        g_total_avl_nodes = 0;
        g_limite_avl = N_INS + 10;
        NoAVL *t = NULL;
        ins_avl = cronoNs(N_INS, [&]
                          { for (int i=0;i<N_INS;i++) t = inserirAVL(t, areas[i], (unsigned)i, "bench"); });
        liberarAVL(t);
    }

    TabelaHash *t_hc = criarTabelaHash();
    double ins_hc = cronoNs(N_INS, [&]
                            { for (int i=0;i<N_INS;i++) inserirHash(t_hc, (unsigned)i, areas[i]); });
    liberarTabelaHash(t_hc);

    TabelaHashOtm *t_ho = criarHashOtm();
    double ins_ho = cronoNs(N_INS, [&]
                            { for (int i=0;i<N_INS;i++) inserirHashOtm(t_ho, (unsigned)i, areas[i]); });
    liberarHashOtm(t_ho);

    NoKD *t_kd = NULL;
    double ins_kd = cronoNs(N_INS, [&]
                            { for (int i=0;i<N_INS;i++){ double pt[2]={areas[i], areas[(i*7)%N_INS]}; t_kd = inserirKD(t_kd, pt, (unsigned)i, 0); } });
    liberarKD(t_kd);

    BloomFilter *t_bl = criarBloom();
    double ins_bl = cronoNs(N_INS, [&]
                            { for (int i=0;i<N_INS;i++) inserirBloom(t_bl, (unsigned)i); });
    liberarBloom(t_bl);

    double mx_i = max({ins_lst, ins_avl, ins_hc, ins_ho, ins_kd, ins_bl, 0.01});
    printf(BOLD "\n  [ INSERCAO - todas as estruturas (%d elementos cada) ]\n" RST, N_INS);
    printLinha('-', 57);
    barraComparativa("Lista (prepend O(1)):", ins_lst, mx_i);
    barraComparativa("AVL O(log n):", ins_avl, mx_i);
    barraComparativa("Hash Encadeamento O(1):", ins_hc, mx_i);
    barraComparativa("Hash Open Addr. O(1):", ins_ho, mx_i);
    barraComparativa("KD-Tree O(log n):", ins_kd, mx_i);
    barraComparativa("Bloom Filter O(k):", ins_bl, mx_i);

    // ===== REMOCAO: estruturas que suportam (Bloom e KD = N/A justificado) =====
    const int N_REM = 10000;
    double rem_lst = 0, rem_avl = 0, rem_hc = 0, rem_ho = 0;
    // IDs pseudo-aleatorios: o compilador nao consegue simular as estruturas em
    vector<unsigned int> ids(N_REM);
    {
        mt19937 g(123);
        for (auto &x : ids)
            x = g();
    }

    {
        NoLista *t = NULL;
        for (int i = 0; i < N_REM; i++)
            t = inserirLista(t, areas[i], ids[i]);
        rem_lst = cronoNs(N_REM, [&]
                          { for (int i=0;i<N_REM;i++) t = removerLista(t, ids[i]); });
        liberarLista(t);
    }

    {
        EstadoAVLGuard g;
        g_total_avl_nodes = 0;
        g_limite_avl = N_REM + 10;
        NoAVL *t = NULL;
        for (int i = 0; i < N_REM; i++)
            t = inserirAVL(t, areas[i], ids[i], "bench");
        rem_avl = cronoNs(N_REM, [&]
                          { for (int i=0;i<N_REM;i++) t = removerAVL(t, areas[i], ids[i]); });
        liberarAVL(t);
    }

    // Hash: medida em buffers estaticos (medirRemChain/medirRemOpen). Usa MUITOS
    const int N_HREM = 300000;
    vector<unsigned int> hids(N_HREM);
    vector<double> hareas(N_HREM);
    {
        mt19937 g(123);
        uniform_real_distribution<double> d(0.001, 0.999);
        for (int i = 0; i < N_HREM; i++)
        {
            hids[i] = g();
            hareas[i] = d(g);
        }
    }
    rem_hc = medirRemChain(hids, hareas);
    rem_ho = medirRemOpen(hids, hareas);

    double mx_r = max({rem_lst, rem_avl, rem_hc, rem_ho, 0.01});
    printf(BOLD "\n  [ REMOCAO - estruturas que suportam ]\n" RST);
    printLinha('-', 57);
    barraComparativa("Lista O(n):", rem_lst, mx_r);
    barraComparativa("AVL O(log n):", rem_avl, mx_r);
    barraComparativa("Hash Encadeamento O(1):", rem_hc, mx_r);
    barraComparativa("Hash Open Addr. O(1):", rem_ho, mx_r);
    printf("  %-28s " YEL "N/A" RST DIM " - Bloom nao remove (apagar bits afetaria outros ids)\n" RST, "Bloom Filter:");
    printf("  %-28s " YEL "N/A" RST DIM " - KD-Tree: remocao exige reconstruir a subarvore\n" RST, "KD-Tree:");

    // MEMÓRIA ESTIMADA
    long long mem_henc = (long long)TAMANHO_TABELA * sizeof(NoHash *) + sizeof(TabelaHash);
    long long mem_hotm = (long long)TAM_HASH_OTM * sizeof(HashSlot);
    long long mem_blm = (TAM_BLOOM + 7) / 8;
    long long mem_avl = (long long)g_total_avl_nodes * sizeof(NoAVL);
    long long mx_mem = max({mem_henc, mem_hotm, mem_blm, mem_avl});

    printf(BOLD "\n  [ USO DE MEMORIA POR ESTRUTURA (estimado via sizeof) ]\n" RST);
    printLinha('-', 57);
    auto barraKB = [&](const char *n, long long bytes)
    {
        printf("  %-28s ", n);
        printBarra((double)bytes, (double)mx_mem, 28, BLU);
        printf(" %6lld KB\n", bytes / 1024);
    };
    barraKB("Hash Encadeamento (tabela):", mem_henc);
    barraKB("Hash Open Addressing:", mem_hotm);
    barraKB("Bloom Filter:", mem_blm);
    barraKB("AVL (nos alocados):", mem_avl);

    // MEMORIA REAL do processo (medida pelo SO, nao estimada)
    long long rss = memoriaRSS_KB();
    printf(BOLD "\n  [ MEMORIA RAM REAL DO PROCESSO (medida pelo SO) ]\n" RST);
    printLinha('-', 57);
    if (rss >= 0)
    {
        printf("  Working Set (RSS) atual: " GRN "%lld KB" RST " (~%.1f MB)\n",
               rss, rss / 1024.0);
        long long soma_est = (mem_henc + mem_hotm + mem_blm + mem_avl) / 1024;
        printf("  Soma estimada (sizeof):  " YEL "%lld KB\n" RST, soma_est);
        printf(DIM "  RSS conta so paginas RESIDENTES. As tabelas hash sao zeradas:\n"
                   "  slots vazios (nunca escritos) nao entram no working set por\n"
                   "  causa da alocacao preguicosa do SO -> RSS pode ser < estimativa.\n" RST);
    }
    else
    {
        printf(YEL "  (medicao de RSS indisponivel nesta plataforma)\n" RST);
    }

    // COLISÕES + fator de carga explicado
    long long n_elem = g_total_avl_nodes;
    double fator_carga = (TAMANHO_TABELA > 0) ? (double)n_elem / TAMANHO_TABELA * 100.0 : 0.0;

    // Demo de colisao com tabela pequena (mostra o fenomeno em escala)
    const int DEMO_TAM = 97;
    long long colisoes_demo = 0;
    {
        vector<NoHash *> demo(DEMO_TAM, nullptr);
        mt19937 g(7);
        uniform_int_distribution<int> dID(1, 500);
        for (int i = 0; i < 200; i++)
        {
            int id = dID(g);
            int slot = id % DEMO_TAM;
            if (demo[slot])
                colisoes_demo++;
            // nao aloca de verdade, so conta
        }
    }

    printf(BOLD "\n  [ COLISOES - Hash Encadeamento ]\n" RST);
    printLinha('-', 57);
    printf("  Tabela principal (%d slots, %lld elementos):\n", TAMANHO_TABELA, n_elem);
    printf("  Fator de carga:     " YEL "%.2f%%\n" RST, fator_carga);
    printf("  Colisoes reais:     " YEL "%lld\n" RST, hash->total_colisoes);
    if (fator_carga < 10.0)
        printf(DIM "  (Fator de carga baixo -> colisoes raras. Normal para este dataset.)\n" RST);
    printf("\n  Demo colisao (tabela de %d slots, 200 insercoes):\n", DEMO_TAM);
    printf("  Colisoes esperadas: " YEL "%lld (fator de carga ~206%%)\n" RST, colisoes_demo);
    printf(DIM "  Quanto maior o fator de carga, mais colisoes e mais lenta a busca.\n" RST);

    // LATÊNCIA MÉDIA DE OPERAÇÕES COMBINADAS (metrica exigida pelo PDF)
    const int N_CICLOS = 50000;
    auto t_ciclo1 = chrono::high_resolution_clock::now();
    for (int i = 0; i < N_CICLOS; i++)
    {
        unsigned int id = 7000000 + i;
        double a = 0.033;
        inserirHash(hash, id, a);
        g_sink_b = buscarHash(hash, id, &a);
        removerHash(hash, id);
    }
    long long lat_henc = chrono::duration_cast<chrono::nanoseconds>(
                             chrono::high_resolution_clock::now() - t_ciclo1)
                             .count() /
                         N_CICLOS;

    auto t_ciclo2 = chrono::high_resolution_clock::now();
    for (int i = 0; i < N_CICLOS; i++)
    {
        unsigned int id = 7000000 + i;
        double a = 0.033;
        inserirHashOtm(hotm, id, a);
        g_sink_b = buscarHashOtm(hotm, id, &a);
        removerHashOtm(hotm, id);
    }
    long long lat_hotm = chrono::duration_cast<chrono::nanoseconds>(
                             chrono::high_resolution_clock::now() - t_ciclo2)
                             .count() /
                         N_CICLOS;

    printf(BOLD "\n  [ LATENCIA MEDIA - operacoes combinadas (ins+busca+rem) ]\n" RST);
    printLinha('-', 57);
    long long mx_lat = max(lat_henc, lat_hotm);
    if (mx_lat == 0)
        mx_lat = 1;
    printf("  %-28s ", "Hash Encadeamento:");
    printBarra(lat_henc, mx_lat, 28, CYN);
    printf(" %lld ns/ciclo\n", lat_henc);
    printf("  %-28s ", "Hash Open Addr. (otim.):");
    printBarra(lat_hotm, mx_lat, 28, CYN);
    printf(" %lld ns/ciclo\n", lat_hotm);

    printLinha('=', 57);
    printf(DIM "  Complexidade teorica: Hash O(1) | AVL O(log n) | Lista O(n)\n" RST);

    // EXPORTA RESULTADOS PARA CSV (para graficos do relatorio)
    ofstream csv("resultados_benchmark.csv");
    if (csv.is_open())
    {
        csv << "operacao,estrutura,ns_por_op\n";
        csv << "busca,hash_encadeamento," << ns_henc << "\n";
        csv << "busca,hash_open_addr," << ns_hotm << "\n";
        csv << "busca,avl," << ns_avl << "\n";
        csv << "busca,lista," << ns_lst << "\n";
        csv << "busca,bloom," << ns_blm << "\n";
        csv << "insercao,lista," << ins_lst << "\n";
        csv << "insercao,avl," << ins_avl << "\n";
        csv << "insercao,hash_encadeamento," << ins_hc << "\n";
        csv << "insercao,hash_open_addr," << ins_ho << "\n";
        csv << "insercao,kd_tree," << ins_kd << "\n";
        csv << "insercao,bloom," << ins_bl << "\n";
        csv << "remocao,lista," << rem_lst << "\n";
        csv << "remocao,avl," << rem_avl << "\n";
        csv << "remocao,hash_encadeamento," << rem_hc << "\n";
        csv << "remocao,hash_open_addr," << rem_ho << "\n";
        csv << "latencia_combinada,hash_encadeamento," << lat_henc << "\n";
        csv << "latencia_combinada,hash_open_addr," << lat_hotm << "\n";
        csv << "memoria_kb,hash_encadeamento," << mem_henc / 1024 << "\n";
        csv << "memoria_kb,hash_open_addr," << mem_hotm / 1024 << "\n";
        csv << "memoria_kb,bloom," << mem_blm / 1024 << "\n";
        csv << "memoria_kb,avl," << mem_avl / 1024 << "\n";
        if (rss >= 0)
            csv << "memoria_kb,rss_processo_real," << rss << "\n";
        csv.close();
        printf(GRN "  [CSV] Resultados salvos em 'resultados_benchmark.csv'\n" RST);
    }
}

// =============================================================================
//  12. TESTE DE ESCALABILIDADE (O(n) vs O(log n) vs O(1))
// =============================================================================

void testeEscalabilidade()
{
    printf("\n" CYN BOLD "  ESCALABILIDADE - Tempo medio de Busca por N\n" RST);
    printLinha('=', 57);
    printf(BOLD "  %-8s  %14s %14s %14s\n" RST,
           "N", "Lista O(n) ns", "AVL O(lgN) ns", "Hash O(1) ns");
    printLinha('-', 57);

    int tamanhos[] = {1000, 5000, 10000, 50000, 100000};
    const int NTAM = 5;
    double ns_lst[NTAM], ns_avl[NTAM], ns_hsh[NTAM];

    for (int t = 0; t < NTAM; t++)
    {
        int N = tamanhos[t];
        NoLista *lst = NULL;
        NoAVL *avl_t = NULL;
        TabelaHash *hsh = criarTabelaHash();
        EstadoAVLGuard guard;

        mt19937 gen(t * 31 + 7);
        uniform_real_distribution<double> d(0.001, 0.999);

        for (int i = 1; i <= N; i++)
        {
            double a = d(gen);
            lst = inserirLista(lst, a, i);
            avl_t = inserirAVL(avl_t, a, i, "esc");
            inserirHash(hsh, i, a);
        }

        // Coleta 32 pares (area, id) distribuidos para variar a busca entre iteracoes.
        const int NP = 32;
        double pa[NP];
        unsigned int pi_arr[NP];
        {
            int stp = max(1, N / NP);
            for (int p = 0; p < NP; p++)
            {
                pi_arr[p] = (unsigned int)(1 + p * stp);
                pa[p] = 0.0;
                buscarHash(hsh, pi_arr[p], &pa[p]);
            }
        }

        // Repeticoes ADAPTATIVAS: a lista e O(n) (cara), AVL/Hash sao baratas.
        const int REP_LENTO = 3000;
        const int REP_RAPIDO = 2000000;
        // volatile counter impede o compilador de hoistar chamadas puras
        auto medirBloco = [&](int reps, auto fn) -> double
        {
            volatile int ctr = 0;
            auto t1 = chrono::high_resolution_clock::now();
            for (int i = 0; i < reps; i++)
            {
                int idx = (ctr++) % NP;
                fn(pa[idx], pi_arr[idx]);
            }
            long long tot = chrono::duration_cast<chrono::nanoseconds>(
                                chrono::high_resolution_clock::now() - t1)
                                .count();
            return (double)tot / reps;
        };

        ns_lst[t] = medirBloco(REP_LENTO, [&](double a, unsigned int)
                               { g_sink_b = buscarLista(lst, a); });
        ns_avl[t] = medirBloco(REP_RAPIDO, [&](double a, unsigned int)
                               { g_sink_b = (buscarAVL(avl_t, a) != NULL); });
        double ar = 0;
        ns_hsh[t] = medirBloco(REP_RAPIDO, [&](double, unsigned int id)
                               { g_sink_b = buscarHash(hsh, id, &ar); g_sink_d = ar; });

        // exibe ns/op com 2 decimais (revela valores sub-nanosegundo sem truncar p/ 0)
        printf("  %-8d  %11.2f ns %11.2f ns %11.2f ns\n",
               N, ns_lst[t], ns_avl[t], ns_hsh[t]);

        liberarLista(lst);
        liberarAVL(avl_t);
        liberarTabelaHash(hsh);
        // guard restaura g_total_avl_nodes automaticamente aqui
    }

    // Grafico proporcional para N=100k
    double mx = max({ns_lst[4], ns_avl[4], ns_hsh[4], 0.01});

    printf(BOLD "\n  Comparacao N=100.000 (proporcional ao mais lento):\n" RST);
    printf("  Lista  ");
    printBarra(ns_lst[4], mx, 34, RED);
    printf(" %.2f ns\n", ns_lst[4]);
    printf("  AVL    ");
    printBarra(max(ns_avl[4], 0.01), mx, 34, YEL);
    printf(" %.2f ns\n", ns_avl[4]);
    printf("  Hash   ");
    printBarra(max(ns_hsh[4], 0.01), mx, 34, GRN);
    printf(" %.2f ns\n", ns_hsh[4]);

    // Razoes de crescimento (usa N=1k e N=100k)
    printf(BOLD "\n  Razao de crescimento N=1k -> N=100k:\n" RST);
    double r_lst = ns_lst[0] > 0 ? ns_lst[4] / ns_lst[0] : 0;
    double r_avl = ns_avl[0] > 0 ? ns_avl[4] / ns_avl[0] : 0;
    double r_hsh = ns_hsh[0] > 0 ? ns_hsh[4] / ns_hsh[0] : 0;
    printf("  Lista: " RED "%.1fx" RST " (esperado ~100x para O(n))\n", r_lst);
    printf("  AVL:   " YEL "%.1fx" RST " (esperado ~%.1fx para O(log n))\n",
           r_avl, log2(100000.0) / log2(1000.0));
    printf("  Hash:  " GRN "%.1fx" RST " (esperado ~1x  para O(1))\n", r_hsh);
    printf(DIM "  Nota: AVL e Hash sub-microsegundo = dados em cache L1. Diferenca\n"
               "  real aparece com N>>1M ou acesso a disco (B-Tree vs lista).\n" RST);

    // EXPORTA para CSV (linha por N -> ideal para grafico de linha no relatorio)
    ofstream csv("resultados_escalabilidade.csv");
    if (csv.is_open())
    {
        csv << "N,lista_ns,avl_ns,hash_ns\n";
        for (int t = 0; t < NTAM; t++)
            csv << tamanhos[t] << "," << ns_lst[t] << "," << ns_avl[t] << "," << ns_hsh[t] << "\n";
        csv.close();
        printf(GRN "  [CSV] Dados salvos em 'resultados_escalabilidade.csv'\n" RST);
    }

    printLinha('=', 57);
    printf(DIM "  Busca no meio da lista = pior caso para O(n), melhor caso para O(log n)/O(1)\n" RST);
}

// =============================================================================
//  12b. DEMONSTRACAO DA RESTRICAO R2 (teto de memoria da AVL)
// =============================================================================

void testarRestricaoR2()
{
    printf("\n" CYN BOLD "  RESTRICAO R2 - Teto de Memoria da AVL\n" RST);
    printLinha('=', 57);
    printf("  O PDF (categoria 1) exige limitar o tamanho de uma estrutura.\n");
    printf("  Aqui a AVL tem um teto rigido: insercoes alem do limite sao\n");
    printf("  rejeitadas para nao ultrapassar o orcamento de memoria.\n");
    printLinha('-', 57);

    // Demonstra em uma AVL ISOLADA (nao mexe na arvore principal)
    const int TETO_DEMO = 5000;
    const int TENTATIVAS = 8000;

    // O guard salva o estado real da AVL e o restaura ao fim da funcao (RAII)
    EstadoAVLGuard guard;
    int teto_real = guard.limite;
    g_total_avl_nodes = 0;
    g_limite_avl = TETO_DEMO;
    g_avl_rejeitados = 0;

    printf("  Teto configurado:     " YEL "%d nos\n" RST, TETO_DEMO);
    printf("  Insercoes tentadas:   " YEL "%d\n" RST, TENTATIVAS);

    NoAVL *demo = NULL;
    mt19937 gen(123);
    uniform_real_distribution<double> d(0.001, 0.999);
    for (int i = 1; i <= TENTATIVAS; i++)
        demo = inserirAVL(demo, d(gen), i, "r2_demo");

    int aceitos = g_total_avl_nodes;
    long long rej = g_avl_rejeitados;

    printLinha('-', 57);
    printf("  Nos aceitos:          " GRN "%d" RST "  ", aceitos);
    printBarra(aceitos, TENTATIVAS, 24, GRN);
    printf("\n");
    printf("  Insercoes rejeitadas: " RED "%lld" RST "  ", rej);
    printBarra((double)rej, TENTATIVAS, 24, RED);
    printf("\n");
    printLinha('-', 57);

    if (aceitos == TETO_DEMO && rej == TENTATIVAS - TETO_DEMO)
        printf(GRN "  [OK] Restricao R2 funcionando: a AVL travou exatamente em %d nos.\n" RST, TETO_DEMO);
    else
        printf(YEL "  [!] Resultado inesperado (verifique a logica do teto).\n" RST);

    liberarAVL(demo);

    printLinha('=', 57);
    printf(DIM "  A arvore principal nao foi afetada (teste isolado). Teto real: %d nos.\n" RST,
           teto_real);
}


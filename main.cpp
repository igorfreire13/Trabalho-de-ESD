#include "comum.h"
#include "util.h"
#include "estruturas.h"
#include "analise.h"
#include "benchmark.h"
#include <cstdio>
#include <iostream>
using namespace std;

// Definicao dos globais declarados em comum.h
int g_total_avl_nodes = 0;
bool g_janela_aberta = true;
int g_limite_avl = LIMITE_MAX_AVL;
long long g_avl_rejeitados = 0;
volatile bool g_sink_b = false;
volatile double g_sink_d = 0.0;

// =============================================================================
//  13. MENU PRINCIPAL
// =============================================================================

int main()
{
    habilitarANSI();

    NoLista *lista = NULL;
    TabelaHashOtm *hotm = criarHashOtm();
    NoAVL *avl = NULL;
    TabelaHash *hash = criarTabelaHash();
    NoKD *kd = NULL;
    BloomFilter *bloom = criarBloom();

    int op = 0;
    do
    {
        printf("\n");
        printHeader();
        printf(GRN "  >> DEMONSTRACAO: execute as opcoes 1 a 8 NESTA ORDEM <<\n" RST);
        printf(DIM "  --- 1) Carregar e tratar os dados -------------------\n" RST);
        printf(WHT "   1." RST " Ingerir Dataset            " DIM "(R11 rede, R18 sensor)\n" RST);
        printf(WHT "   2." RST " Histograma de Distribuicao " DIM "(ver dados brutos)\n" RST);
        printf(WHT "   3." RST " Filtrar Outliers Z/IQR     " DIM "(R2 limite, R18 poda)\n" RST);
        printf(DIM "  --- 2) Analisar -------------------------------------\n" RST);
        printf(WHT "   4." RST " Analise Espacial KD-Tree   " DIM "(produtividade por frame)\n" RST);
        printf(DIM "  --- 3) Medir desempenho -----------------------------\n" RST);
        printf(WHT "   5." RST " Teste de Escalabilidade    " DIM "(+CSV p/ relatorio)\n" RST);
        printf(WHT "   6." RST " Benchmark Completo         " DIM "(+CSV p/ relatorio)\n" RST);
        printf(DIM "  --- 4) Demonstrar restricoes ------------------------\n" RST);
        printf(WHT "   7." RST " Teste do Teto da AVL       " DIM "(restricao R2 memoria)\n" RST);
        printf(WHT "   8." RST " Teste de Bloqueio de CPU   " DIM "(restricao R9)\n" RST);
        printf(DIM "  --- Operacoes avulsas (opcional) --------------------\n" RST);
        printf(WHT "   9." RST " Buscar registro por ID\n");
        printf(WHT "  10." RST " Inserir registro individual\n");
        printf(WHT "  11." RST " Remover registro por ID\n");
        printf(DIM "  -----------------------------------------------------\n" RST);
        printf(WHT "   0." RST " Encerrar e liberar memoria\n");
        printLinha('-');
        printf(CYN "  Escolha: " RST);

        if (!(cin >> op))
        {
            cin.clear();
            cin.ignore(10000, '\n');
            op = -1;
        }

        switch (op)
        {

        case 1:
            carregarDataset(&lista, hotm, &avl, hash, &kd, bloom);
            printf(YEL "  Status AVL (R2): " RST "%d / %d nos (%.1f%% do teto)",
                   g_total_avl_nodes, g_limite_avl,
                   100.0 * g_total_avl_nodes / g_limite_avl);
            if (g_avl_rejeitados > 0)
                printf(RED "  | %lld rejeitadas pelo teto!\n" RST, g_avl_rejeitados);
            else
                printf(DIM "  | teto nao atingido (use a opcao 7 p/ ver R2 ativo)\n" RST);
            break;

        case 4:
            if (!kd)
            {
                printf(RED "  [!] Ingira os dados primeiro (opcao 1).\n" RST);
                break;
            }
            analiseEspacialKD(kd);
            break;

        case 3:
        {
            if (!lista)
            {
                printf(RED "  [!] Ingira os dados primeiro (opcao 1).\n" RST);
                break;
            }

            vector<double> vals;
            for (NoLista *p = lista; p; p = p->prox)
                vals.push_back(p->bounding_area);
            ResumoEstatistico r = calcularResumo(vals);

            // Conta quantas anomalias cada metodo capturaria
            auto contarFora = [&](double lo, double hi)
            {
                int c = 0;
                for (double v : vals)
                    if (v < lo || v > hi)
                        c++;
                return c;
            };
            int n_zsc = contarFora(r.zsc_inf, r.zsc_sup);
            int n_iqr = contarFora(r.iqr_inf, r.iqr_sup);

            printf("\n" CYN BOLD "  DETECCAO DE OUTLIERS - Z-SCORE vs IQR\n" RST);
            printLinha('-', 57);
            printf("  Z-Score (media +/- 3 sigma):\n");
            printf("    Faixa: " YEL "[%.4f , %.4f]" RST "  capturaria " RED "%d" RST " outliers\n",
                   r.zsc_inf, r.zsc_sup, n_zsc);
            printf("    " DIM "Problema: as anomalias inflam sigma (%.4f), alargando a faixa.\n" RST, r.desvio);
            printf("  IQR / Tukey (Q3 + 1.5*IQR):\n");
            printf("    Faixa: " YEL "[%.4f , %.4f]" RST "  capturaria " GRN "%d" RST " outliers\n",
                   r.iqr_inf, r.iqr_sup, n_iqr);
            printf("    " DIM "Robusto: baseado em mediana/quartis, imune a contaminacao.\n" RST);
            printLinha('-', 57);

            // Aplica o metodo robusto (IQR) na AVL
            printf("  Aplicando filtro IQR (robusto) na AVL (R18)...\n");
            int antes = g_total_avl_nodes;
            avl = podarAVL(avl, r.iqr_inf, r.iqr_sup);
            printf(GRN "  [OK] %d outliers removidos. Nos validos na AVL: %d\n" RST,
                   antes - g_total_avl_nodes, g_total_avl_nodes);
            printf(DIM "  Compare: o IQR removeu mais anomalias R18 que o z-score teria removido.\n" RST);
            printLinha('-', 57);
            break;
        }

        case 6:
            if (!lista)
            {
                printf(RED "  [!] Ingira os dados primeiro (opcao 1).\n" RST);
                break;
            }
            executarBenchmark(lista, hotm, avl, hash, bloom, lista->sample_id);
            break;

        case 5:
            testeEscalabilidade();
            break;

        case 2:
            exibirHistograma(lista, avl);
            break;

        case 8:
        {
            double ar = 0.0;
            unsigned int id_t = lista ? lista->sample_id : 1050;
            printf("\n" CYN BOLD "  TESTE R9 -Disponibilidade da CPU\n" RST);
            printLinha('-', 55);
            g_janela_aberta = false;
            printf("  Passo 1: " RED "CPU 100%% ocupada" RST " -tentando buscar ID %u:\n", id_t);
            buscarHashRestrita(hash, id_t, &ar);
            g_janela_aberta = true;
            printf("\n  Passo 2: " GRN "CPU liberada" RST " -tentando buscar ID %u:\n", id_t);
            if (buscarHashRestrita(hash, id_t, &ar))
                printf(GRN "  [OK] ID %u encontrado. Area: %.6f\n" RST, id_t, ar);
            else
                printf(YEL "  [!] ID %u nao encontrado (sem dados carregados).\n" RST, id_t);
            printLinha('-', 55);
            break;
        }

        case 10:
        {
            unsigned int nid = (unsigned int)lerInteiro("  ID do novo registro: ");
            double narea = lerDouble("  Bounding area:       ");
            lista = inserirLista(lista, narea, nid);
            inserirHashOtm(hotm, nid, narea);
            avl = inserirAVL(avl, narea, nid, "manual.png");
            inserirHash(hash, nid, narea);
            inserirBloom(bloom, nid);
            printf(GRN "  [OK] ID %u inserido em todas as estruturas.\n" RST, nid);
            break;
        }

        case 9:
        {
            unsigned int bid = (unsigned int)lerInteiro("  ID para buscar: ");
            double bar = 0.0;
            bool bBlm = buscarBloom(bloom, bid);
            bool bHsh = buscarHash(hash, bid, &bar);
            double bOptArea = 0.0;
            bool bOpt = buscarHashOtm(hotm, bid, &bOptArea);
            NoAVL *bAvl = buscarAVL(avl, bar);

            printf("\n" CYN BOLD "  RESULTADO -ID %u\n" RST, bid);
            printLinha('-', 55);
            // Bloom nao suporta delecao: pode retornar true mesmo apos remocao (falso positivo)
            if (bBlm && !bHsh)
                printf("  Bloom Filter:    " YEL "FALSO POSITIVO (elemento foi removido, Bloom nao suporta delecao)\n" RST);
            else if (bBlm)
                printf("  Bloom Filter:    " GRN "Possivel membro (confirmar nas estruturas exatas abaixo)\n" RST);
            else
                printf("  Bloom Filter:    " RED "AUSENTE com certeza (nenhum falso negativo possivel)\n" RST);
            printf("  Hash Encadeada:  %s  area=%.6f\n", bHsh ? GRN "ENCONTRADO" RST : YEL "nao encontrado" RST, bar);
            printf("  Hash Otimizada:  %s  area=%.6f\n", bOpt ? GRN "ENCONTRADO" RST : YEL "nao encontrado" RST, bOptArea);
            printf("  AVL:             %s\n", bAvl ? GRN "ENCONTRADO" RST : YEL "nao encontrado (pode ter sido podado)" RST);
            printLinha('-', 55);
            break;
        }

        case 11:
        {
            unsigned int rid = (unsigned int)lerInteiro("  ID para remover: ");
            double rar = 0.0;
            if (!buscarHash(hash, rid, &rar))
            {
                printf(YEL "  [!] ID %u nao encontrado.\n" RST, rid);
                break;
            }
            lista = removerLista(lista, rid);
            removerHashOtm(hotm, rid);
            avl = removerAVL(avl, rar, rid);
            removerHash(hash, rid);
            printf(GRN "  [OK] ID %u removido de todas as estruturas.\n" RST, rid);
            break;
        }

        case 7:
            testarRestricaoR2();
            break;

        case 0:
            printf(CYN "  Desalocando memoria e encerrando...\n" RST);
            break;

        default:
            printf(RED "  [!] Opcao invalida.\n" RST);
        }
    } while (op != 0);

    liberarLista(lista);
    liberarHashOtm(hotm);
    liberarAVL(avl);
    liberarTabelaHash(hash);
    liberarKD(kd);
    liberarBloom(bloom);
    return 0;
}

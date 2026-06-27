#include "analise.h"
#include "estruturas.h"
#include "util.h"
#include "comum.h"
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <random>
#include <thread>
#include <chrono>
#include <unordered_map>
using namespace std;

// =============================================================================
//  8. OPERACOES ADICIONAIS (3 exigidas pelo PDF)
//  (a) Estatistica (media/desvio/outliers Z-Score vs IQR); (b) filtragem (poda
//  da AVL); (c) classificacao espacial (densidade por frame na KD-Tree).
// =============================================================================

double calcularMedia(NoLista *cab, int *total)
{
    double soma = 0.0;
    int n = 0;
    for (NoLista *p = cab; p; p = p->prox)
    {
        soma += p->bounding_area;
        n++;
    }
    *total = n;
    return n ? soma / n : 0.0;
}
double calcularDP(NoLista *cab, double media, int total)
{
    if (!cab || total <= 1)
        return 0.0;
    double sq = 0.0;
    for (NoLista *p = cab; p; p = p->prox)
    {
        double d = p->bounding_area - media;
        sq += d * d;
    }
    return sqrt(sq / (total - 1));
}

// --- Estatistica robusta: percentil de um vetor JA ORDENADO (interpolacao linear) ---
static double percentil(const vector<double> &ord, double p)
{
    if (ord.empty())
        return 0.0;
    if (ord.size() == 1)
        return ord[0];
    double idx = p / 100.0 * (ord.size() - 1);
    int lo = (int)idx;
    double frac = idx - lo;
    if (lo + 1 >= (int)ord.size())
        return ord.back();
    return ord[lo] * (1 - frac) + ord[lo + 1] * frac;
}
ResumoEstatistico calcularResumo(vector<double> vals)
{
    ResumoEstatistico r = {};
    if (vals.empty())
        return r;
    sort(vals.begin(), vals.end());

    r.minimo = vals.front();
    r.maximo = vals.back();
    r.q1 = percentil(vals, 25);
    r.mediana = percentil(vals, 50);
    r.q3 = percentil(vals, 75);
    r.iqr = r.q3 - r.q1;

    double soma = 0.0;
    for (double v : vals)
        soma += v;
    r.media = soma / vals.size();

    double sq = 0.0;
    for (double v : vals)
    {
        double d = v - r.media;
        sq += d * d;
    }
    r.desvio = vals.size() > 1 ? sqrt(sq / (vals.size() - 1)) : 0.0;

    // Tukey fences (robusto): Q1 - 1.5*IQR  ate  Q3 + 1.5*IQR
    r.iqr_inf = max(0.0, r.q1 - 1.5 * r.iqr);
    r.iqr_sup = r.q3 + 1.5 * r.iqr;

    // Z-score (sensivel a contaminacao): media +/- 3*sigma
    r.zsc_inf = max(0.0, r.media - 3 * r.desvio);
    r.zsc_sup = r.media + 3 * r.desvio;
    return r;
}

static void painelHistograma(const vector<double> &vals, const char *titulo,
                             const char *cor, bool log_escala)
{
    if (vals.empty())
    {
        printf(YEL "  (sem dados)\n" RST);
        return;
    }
    const int NB = 10;
    double mn = *min_element(vals.begin(), vals.end());
    double mx = *max_element(vals.begin(), vals.end());
    double step = (mx - mn);
    if (step < 1e-9)
    {
        printf(YEL "  Todos os valores sao identicos (%.4f)\n" RST, mn);
        return;
    }
    step /= NB;

    int cnt[10] = {};
    for (double v : vals)
    {
        int b = (int)((v - mn) / step);
        if (b >= NB)
            b = NB - 1;
        cnt[b]++;
    }
    int mx_cnt = *max_element(cnt, cnt + NB);

    printf(CYN BOLD "  %s  (%zu valores)%s\n" RST, titulo, vals.size(),
           log_escala ? "  [escala LOG]" : "  [escala linear]");
    printLinha('-', 57);
    for (int i = 0; i < NB; i++)
    {
        double lo = mn + i * step, hi = lo + step;
        printf(YEL "  [%6.4f-%6.4f] " RST, lo, hi);
        int barras;
        if (log_escala)
            // log1p comprime: balde de 1 ainda aparece, balde de 23k nao estoura
            barras = (cnt[i] > 0) ? (int)(log1p(cnt[i]) / log1p(mx_cnt) * 34) : 0;
        else
            barras = mx_cnt > 0 ? cnt[i] * 34 / mx_cnt : 0;
        printf("%s", cor);
        for (int j = 0; j < barras; j++)
            putchar('#');
        printf(RST " %d\n", cnt[i]);
    }
    printLinha('-', 57);
}

// --- Painel de resumo estatistico (explica a assimetria) ---
static void painelResumo(const ResumoEstatistico &r)
{
    printf(CYN BOLD "  RESUMO ESTATISTICO\n" RST);
    printLinha('-', 57);
    printf("  Minimo : " YEL "%8.4f" RST "    Q1      : " YEL "%8.4f\n" RST, r.minimo, r.q1);
    printf("  Mediana: " YEL "%8.4f" RST "    Media   : " YEL "%8.4f\n" RST, r.mediana, r.media);
    printf("  Q3     : " YEL "%8.4f" RST "    Maximo  : " YEL "%8.4f\n" RST, r.q3, r.maximo);
    printf("  IQR    : " YEL "%8.4f" RST "    Desvio  : " YEL "%8.4f\n" RST, r.iqr, r.desvio);
    printLinha('-', 57);
    // Diagnostico de assimetria: media vs mediana
    if (r.media > r.mediana * 1.5)
        printf(DIM "  Media (%.4f) >> Mediana (%.4f): distribuicao ASSIMETRICA A DIREITA.\n"
                   "  Causa: ~10%% de anomalias R18 puxam a media para cima, mas nao a mediana.\n" RST,
               r.media, r.mediana);
    printLinha('-', 57);
}

void exibirHistograma(NoLista *cab, NoAVL *avl)
{
    if (!cab)
    {
        printf(RED "  Sem dados. Ingira o dataset primeiro.\n" RST);
        return;
    }

    vector<double> vals_lista;
    for (NoLista *p = cab; p; p = p->prox)
        vals_lista.push_back(p->bounding_area);

    // Resumo estatistico primeiro: explica a assimetria
    printf("\n");
    ResumoEstatistico r = calcularResumo(vals_lista);
    painelResumo(r);

    // Painel 1: linear (mostra o pico dominante = dados reais)
    printf("\n");
    painelHistograma(vals_lista, "DADOS BRUTOS - linear", GRN, false);

    // Painel 2: LOG (revela a estrutura escondida = cauda de anomalias R18)
    printf("\n");
    painelHistograma(vals_lista, "DADOS BRUTOS - logaritmico", CYN, true);

    // Painel 3: AVL atual (reflete poda se ja executada)
    vector<double> vals_avl;
    coletarAVL(avl, vals_avl);
    printf("\n");
    if (vals_avl.size() == vals_lista.size())
        painelHistograma(vals_avl, "AVL (poda nao executada ainda)", YEL, true);
    else
        painelHistograma(vals_avl, "AVL (pos-poda - outliers removidos)", BLU, true);

    printf(DIM "  O painel LOG revela que os baldes 'vazios' tem dezenas a centenas\n"
               "  de anomalias R18. Rode a opcao 3 para filtra-las e compare a AVL.\n" RST);
}
// =============================================================================
//  ANALISE ESPACIAL VIA KD-TREE
//  cx,cy sao a posicao da deteccao no quadro da imagem (por frame), nao a
//  geografia do campo. (A) distribuicao, (B) vizinho, (C) produtividade/cluster.
// =============================================================================

// --- Range query: conta pontos dentro de um retangulo, podando subarvores ---
void contarZona(NoKD *r, double x0, double x1, double y0, double y1, int d, int *c)
{
    if (!r)
        return;
    double x = r->ponto[0], y = r->ponto[1];
    if (x >= x0 && x < x1 && y >= y0 && y < y1)
        (*c)++;
    int cd = d % K;
    if (cd == 0)
    {
        if (x0 < x)
            contarZona(r->esq, x0, x1, y0, y1, d + 1, c);
        if (x1 >= x)
            contarZona(r->dir, x0, x1, y0, y1, d + 1, c);
    }
    else
    {
        if (y0 < y)
            contarZona(r->esq, x0, x1, y0, y1, d + 1, c);
        if (y1 >= y)
            contarZona(r->dir, x0, x1, y0, y1, d + 1, c);
    }
}

// (A) Distribuicao espacial das deteccoes NO QUADRO (nao no campo)
void analisarDistribuicaoEspacial(NoKD *kd)
{
    int q[4] = {};
    contarZona(kd, 0.0, 0.5, 0.5, 1.0, 0, &q[0]);
    contarZona(kd, 0.5, 1.0, 0.5, 1.0, 0, &q[1]);
    contarZona(kd, 0.0, 0.5, 0.0, 0.5, 0, &q[2]);
    contarZona(kd, 0.5, 1.0, 0.0, 0.5, 0, &q[3]);

    int total = q[0] + q[1] + q[2] + q[3];
    if (!total)
    {
        printf(YEL "  Nenhum ponto em [0,1]. Verifique coordenadas.\n" RST);
        return;
    }
    double med = total / 4.0;
    const char *nomes[4] = {"Quadro Sup-Esquerda", "Quadro Sup-Direita",
                            "Quadro Inf-Esquerda", "Quadro Inf-Direita"};

    printf("\n" CYN BOLD "  (A) DISTRIBUICAO DAS DETECCOES NO QUADRO (KD-Tree range query)\n" RST);
    printLinha('-', 57);
    int mx = *max_element(q, q + 4);
    for (int i = 0; i < 4; i++)
    {
        const char *cor = (q[i] >= med * 1.2) ? GRN : (q[i] <= med * 0.8) ? RED
                                                                          : YEL;
        printf("  %-20s ", nomes[i]);
        printBarra(q[i], mx, 18, cor);
        printf(" %5d", q[i]);
        if (q[i] >= med * 1.2)
            printf(GRN " [concentracao ALTA]\n" RST);
        else if (q[i] <= med * 0.8)
            printf(RED " [concentracao BAIXA]\n" RST);
        else
            printf(YEL " [normal]\n" RST);
    }
    printLinha('-', 57);
    printf(DIM "  Mede a regiao DO QUADRO onde as uvas aparecem (vies de camera),\n"
               "  nao a geografia do campo - cx,cy sao coordenadas da imagem.\n" RST);
}

// --- (B) Busca do VIZINHO MAIS PROXIMO (nearest-neighbor) na KD-Tree ---------
static double distSqKD(const double a[], const double b[])
{
    double dx = a[0] - b[0], dy = a[1] - b[1];
    return dx * dx + dy * dy;
}
void nnKD(NoKD *r, const double alvo[], int depth,
          double *melhorDist, unsigned int *melhorId, unsigned int exclId)
{
    if (!r)
        return;
    if (r->sample_id != exclId)
    {
        double dd = distSqKD(r->ponto, alvo);
        if (dd < *melhorDist)
        {
            *melhorDist = dd;
            *melhorId = r->sample_id;
        }
    }
    int cd = depth % K;
    double diff = alvo[cd] - r->ponto[cd];
    NoKD *perto = diff < 0 ? r->esq : r->dir;
    NoKD *longe = diff < 0 ? r->dir : r->esq;
    nnKD(perto, alvo, depth + 1, melhorDist, melhorId, exclId);
    // so explora o outro lado se ainda houver chance de algo mais perto la
    if (diff * diff < *melhorDist)
        nnKD(longe, alvo, depth + 1, melhorDist, melhorId, exclId);
}
static void coletarPontosKD(NoKD *r, vector<PontoKD> &v)
{
    if (!r)
        return;
    coletarPontosKD(r->esq, v);
    v.push_back({r->ponto[0], r->ponto[1], r->sample_id, r->frame_id});
    coletarPontosKD(r->dir, v);
}

// (B) Vizinho mais proximo: dado um ponto, acha a deteccao mais proxima O(log n).
void analisarVizinhanca(NoKD *kd)
{
    printf("\n" CYN BOLD "  (B) BUSCA DO VIZINHO MAIS PROXIMO (KD-Tree NN)\n" RST);
    printLinha('-', 57);
    printf("  Dado um ponto (x,y) no quadro, acha a deteccao mais proxima:\n");

    struct Q
    {
        double x, y;
        const char *nome;
    };
    Q consultas[] = {
        {0.50, 0.50, "Centro"},
        {0.10, 0.10, "Sup-esquerda"},
        {0.90, 0.90, "Inf-direita"},
        {0.90, 0.10, "Sup-direita"},
    };
    for (const Q &q : consultas)
    {
        double pt[2] = {q.x, q.y};
        double bd = 1e18;
        unsigned int bid = 0;
        nnKD(kd, pt, 0, &bd, &bid, 0);
        if (bid)
            printf("    (%.2f, %.2f) %-13s -> " YEL "ID %u" RST " a dist %.4f\n",
                   q.x, q.y, q.nome, bid, sqrt(bd));
    }
    printLinha('-', 57);
    printf(DIM "  A poda da KD-Tree torna a busca O(log n) medio, contra O(n) de\n"
               "  comparar com todas as deteccoes uma a uma.\n" RST);
}

// Range query por raio: coleta ids dentro de um circulo (filtrando por frame).
static void rangeKD(NoKD *r, const double alvo[], double eps, int frame, int depth,
                    vector<unsigned int> &out)
{
    if (!r)
        return;
    if ((frame < 0 || r->frame_id == frame) && distSqKD(r->ponto, alvo) <= eps * eps)
        out.push_back(r->sample_id);
    int cd = depth % K;
    double diff = alvo[cd] - r->ponto[cd];
    NoKD *perto = diff < 0 ? r->esq : r->dir;
    NoKD *longe = diff < 0 ? r->dir : r->esq;
    rangeKD(perto, alvo, eps, frame, depth + 1, out);
    if (diff * diff <= eps * eps)
        rangeKD(longe, alvo, eps, frame, depth + 1, out);
}

// (C) Produtividade por frame: conta cachos por frame e classifica na faixa
// agronomica (poucos = subproducao, muitos = sobrecarga).
void produtividadePorFrame(const vector<PontoKD> &pts)
{
    const int IDEAL_MIN = 9;
    const int IDEAL_MAX = 20;

    unordered_map<int, int> cachosPorFrame;
    for (const PontoKD &p : pts)
        cachosPorFrame[p.frame]++;

    int sub = 0, ideal = 0, sobre = 0, maxc = 0;
    double soma = 0;
    for (auto &kv : cachosPorFrame)
    {
        int c = kv.second;
        if (c < IDEAL_MIN)
            sub++;
        else if (c <= IDEAL_MAX)
            ideal++;
        else
            sobre++;
        soma += c;
        maxc = max(maxc, c);
    }
    int nframes = (int)cachosPorFrame.size();
    if (!nframes)
        return;

    printf("\n" CYN BOLD "  (C) PRODUTIVIDADE POR FRAME (cachos por imagem/planta)\n" RST);
    printLinha('-', 57);
    printf("  Frames analisados:  " YEL "%d" RST "   | media: " YEL "%.1f" RST
           " cachos/frame  | max: " YEL "%d\n" RST,
           nframes, soma / nframes, maxc);
    printLinha('-', 57);
    printf("  Classificacao agronomica (faixa ideal %d-%d cachos/frame):\n", IDEAL_MIN, IDEAL_MAX);
    int mx = max({sub, ideal, sobre, 1});
    printf("    " RED "Subproducao (< %2d): " RST, IDEAL_MIN);
    printBarra(sub, mx, 22, RED);
    printf(" %d frames\n", sub);
    printf("    " GRN "Ideal (%2d-%2d):      " RST, IDEAL_MIN, IDEAL_MAX);
    printBarra(ideal, mx, 22, GRN);
    printf(" %d frames\n", ideal);
    printf("    " YEL "Sobrecarga (> %2d):  " RST, IDEAL_MAX);
    printBarra(sobre, mx, 22, YEL);
    printf(" %d frames -> raleio\n", sobre);
    printLinha('-', 57);
    printf(DIM "  Agora cada frame tem coordenadas proprias (frame_id recuperado dos\n"
               "  nomes dos arquivos do Roboflow), entao cachos/frame e uma metrica\n"
               "  agronomica VALIDA - nao mistura imagens diferentes.\n" RST);
}

// (D) DBSCAN por frame: agrupa cachos adjacentes do mesmo frame (KD-Tree faz as
// vizinhancas em O(log n) -> ~O(n log n), contra O(n^2) da forca bruta).
void clusterizarPorFrame(NoKD *kd, const vector<PontoKD> &pts,
                         const unordered_map<unsigned int, int> &idx)
{
    const int MINPTS = 3;
    const double EPS = 0.10;

    int m = (int)pts.size();
    vector<int> label(m, -1);
    int cid = 0;
    for (int i = 0; i < m; i++)
    {
        if (label[i] != -1)
            continue;
        double alvo[2] = {pts[i].x, pts[i].y};
        vector<unsigned int> viz;
        rangeKD(kd, alvo, EPS, pts[i].frame, 0, viz);
        if ((int)viz.size() < MINPTS)
        {
            label[i] = 0;
            continue;
        }
        cid++;
        label[i] = cid;
        vector<int> fila;
        for (unsigned int id : viz)
        {
            int j = idx.at(id);
            if (j != i)
                fila.push_back(j);
        }
        for (size_t s = 0; s < fila.size(); s++)
        {
            int j = fila[s];
            if (label[j] == 0)
                label[j] = cid;
            if (label[j] != -1)
                continue;
            label[j] = cid;
            double a2[2] = {pts[j].x, pts[j].y};
            vector<unsigned int> viz2;
            rangeKD(kd, a2, EPS, pts[j].frame, 0, viz2);
            if ((int)viz2.size() >= MINPTS)
                for (unsigned int id2 : viz2)
                    fila.push_back(idx.at(id2));
        }
    }
    int isolados = 0;
    for (int i = 0; i < m; i++)
        if (label[i] <= 0)
            isolados++;

    printf("\n" CYN BOLD "  (D) GRUPOS DE CACHOS COLADOS POR FRAME (DBSCAN c/ KD-Tree)\n" RST);
    printLinha('-', 57);
    printf("  eps=%.2f minPts=%d (so conecta cachos do MESMO frame)\n", EPS, MINPTS);
    printf("  Grupos de cachos adjacentes: " YEL "%d" RST "   | cachos isolados: " YEL "%d\n" RST,
           cid, isolados);
    printLinha('-', 57);
    printf(DIM "  Filtrando por frame, o DBSCAN agrupa so cachos realmente vizinhos\n"
               "  na MESMA imagem (cachos colados sao reais, nao erro). A KD-Tree faz\n"
               "  cada consulta de vizinhanca em O(log n) -> ~O(n log n) total.\n" RST);
}

// Relatorio espacial completo (A + B + C + D)
void analiseEspacialKD(NoKD *kd)
{
    analisarDistribuicaoEspacial(kd);
    analisarVizinhanca(kd);

    vector<PontoKD> pts;
    coletarPontosKD(kd, pts);
    if (pts.empty())
        return;
    unordered_map<unsigned int, int> idx;
    idx.reserve(pts.size() * 2);
    for (int i = 0; i < (int)pts.size(); i++)
        idx[pts[i].id] = i;

    produtividadePorFrame(pts);
    clusterizarPorFrame(kd, pts, idx);
}

// =============================================================================
//  9. GERADOR DE DATASET SINTÉTICO
// =============================================================================

void gerarDataset(const string &nome, int n)
{
    printf("  Gerando %d amostras sinteticas em '%s'...\n", n, nome.c_str());
    ofstream f(nome);
    if (!f.is_open())
    {
        printf(RED "  [ERRO] Falha ao criar arquivo.\n" RST);
        return;
    }
    // Schema novo: sample_id, frame_id, cx, cy, width, height, area
    f << "sample_id,frame_id,cx,cy,width,height,area\n";
    mt19937 gen(2026);
    uniform_real_distribution<double> dxy(0.05, 0.95), dwh(0.01, 0.15);
    uniform_int_distribution<int> dCachos(5, 20);
    normal_distribution<double> ruido(0.0, 0.003);
    int sid = 0, frame = 0;
    while (sid < n)
    {
        int cachos = dCachos(gen);
        for (int c = 0; c < cachos && sid < n; c++)
        {
            double cx = dxy(gen), cy = dxy(gen), w = dwh(gen), h = dwh(gen);
            double area = max(0.0001, w * h + ruido(gen));
            f << sid << "," << frame << "," << cx << "," << cy << "," << w << "," << h << "," << area << "\n";
            sid++;
        }
        frame++;
    }
    f.close();
    printf(GRN "  [OK] Dataset gerado com sucesso.\n" RST);
}

// =============================================================================
//  10. INGESTÃO DO CSV
// =============================================================================

void carregarDataset(NoLista **lista, TabelaHashOtm *hotm, NoAVL **avl,
                     TabelaHash *hash, NoKD **kd, BloomFilter *bloom)
{
    string nome = "dataset_uvas_consolidado.csv";
    ifstream arq(nome);
    if (!arq.is_open())
    {
        printf(YEL "\n  [AVISO] '%s' nao encontrado.\n" RST, nome.c_str());
        printf("  Gerar dataset sintetico com 500.000 amostras? (s/n): ");
        char r;
        cin >> r;
        if (r == 's' || r == 'S')
        {
            gerarDataset(nome, 500000);
            arq.open(nome);
        }
        if (!arq.is_open())
            return;
    }

    // conta linhas para barra de progresso
    unsigned int total_linhas = 0;
    {
        string ln;
        getline(arq, ln);
        while (getline(arq, ln))
            total_linhas++;
    }
    arq.clear();
    arq.seekg(0);
    string cabecalho;
    getline(arq, cabecalho);

    printf("\n" CYN "  Ingerindo '%s' (%u registros)...\n" RST, nome.c_str(), total_linhas);

    mt19937 gen(42);
    uniform_real_distribution<double> dNoise(0.40, 0.95);
    uniform_int_distribution<int> dDefeito(1, 100);

    unsigned int cnt = 0;
    string linha;
    while (getline(arq, linha))
    {
        if (!linha.empty() && linha.back() == '\r')
            linha.pop_back();
        if (linha.empty())
            continue;
        stringstream ss(linha);
        string cel;
        vector<string> cols;
        while (getline(ss, cel, ','))
            cols.push_back(cel);
        // Schema novo: sample_id, frame_id, cx, cy, width, height, area (7 colunas)
        if (cols.size() < 7)
            continue;
        try
        {
            unsigned int id = stoul(cols[0]);
            int frame = stoi(cols[1]);
            double cx = stod(cols[2]), cy = stod(cols[3]), area = stod(cols[6]);
            if (dDefeito(gen) <= 10)
                area = dNoise(gen);
            if (cnt % 1000 == 0)
                this_thread::sleep_for(chrono::milliseconds(1));

            *lista = inserirLista(*lista, area, id);
            inserirHashOtm(hotm, id, area);
            *avl = inserirAVL(*avl, area, id, "frame_real.png");
            inserirHash(hash, id, area);
            double pt[2] = {cx, cy};
            *kd = inserirKD(*kd, pt, id, frame);
            inserirBloom(bloom, id);
            cnt++;
            if (cnt % 10000 == 0)
                printProgresso(cnt, total_linhas);
        }
        catch (...)
        {
            continue;
        }
    }
    arq.close();
    printProgresso(cnt, total_linhas);
    printf("\n" GRN "\n  [OK] %u amostras carregadas.\n" RST, cnt);
}

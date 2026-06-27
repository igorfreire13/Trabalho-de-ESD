# Sistema de Análise de Detecções de Uva

Projeto de Estruturas de Dados 2026 — implementação e comparação de estruturas
de dados aplicadas a um dataset de detecções de cachos de uva (visão computacional).

## Estruturas implementadas

- **Lista encadeada** — baseline O(n)
- **Árvore AVL** — busca ordenada O(log n) auto-balanceada
- **Tabela Hash (encadeamento)** — acesso O(1) médio
- **Tabela Hash (open addressing)** — versão otimizada (requisito 7)
- **KD-Tree** — busca espacial 2D por frame (fora da ementa)
- **Bloom Filter** — teste de pertinência O(k) (fora da ementa)

## Organização do código

| Arquivo | Responsabilidade |
|---------|------------------|
| `comum.h` | Constantes, macros, structs e globais compartilhados |
| `util.h/.cpp` | Utilidades de terminal (cores, barras, entrada) |
| `estruturas.h/.cpp` | As 6 estruturas de dados |
| `analise.h/.cpp` | Estatística, histograma, análise espacial e ingestão |
| `benchmark.h/.cpp` | Benchmarks, escalabilidade e teste de restrição |
| `main.cpp` | Menu principal |

## Como compilar e executar

```bash
g++ -O2 -std=c++17 -o sistema util.cpp estruturas.cpp analise.cpp benchmark.cpp main.cpp
./sistema        # no Windows: sistema.exe
```

## Uso

No menu, execute as opções **1 a 8 nesta ordem** para a demonstração completa.
As opções 9 a 11 são operações individuais (buscar/inserir/remover).

## Dados

- `dataset_uvas_consolidado.csv` — dataset usado pelo programa. Foi gerado a
  partir das detecções YOLO (Roboflow), consolidando os arquivos de label em um
  único CSV com `frame_id`. O processo está descrito no relatório.
- `resultados_benchmark.csv` / `resultados_escalabilidade.csv` — saídas dos
  benchmarks (geradas pelo programa, usadas nos gráficos do relatório).

## Restrições implementadas (5 categorias)

R2 (memória), R9 (processamento), R11 (latência), R18 (dados), R21 (algorítmica).

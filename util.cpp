#include "util.h"
#include "comum.h"
#include <cstdio>
#include <iostream>
#ifdef _WIN32
#include <windows.h>
#include <psapi.h>
extern "C" BOOL WINAPI K32GetProcessMemoryInfo(HANDLE, PROCESS_MEMORY_COUNTERS *, DWORD);
#endif
using namespace std;

void habilitarANSI()
{
#ifdef _WIN32
    HANDLE h = GetStdHandle(STD_OUTPUT_HANDLE);
    DWORD m = 0;
    if (GetConsoleMode(h, &m))
        SetConsoleMode(h, m | ENABLE_VIRTUAL_TERMINAL_PROCESSING);
    SetConsoleOutputCP(CP_UTF8);
#endif
}

void printLinha(char c, int n)
{
    printf(DIM);
    for (int i = 0; i < n; i++)
        putchar(c);
    printf(RST "\n");
}

// Retorna a memoria RAM (Working Set / RSS) realmente usada pelo processo, em KB.
long long memoriaRSS_KB()
{
#ifdef _WIN32
    PROCESS_MEMORY_COUNTERS pmc;
    if (K32GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc)))
        return (long long)(pmc.WorkingSetSize / 1024);
#endif
    return -1;
}

// --- Leitura robusta de input: nunca trava o menu com entrada invalida ---
int lerInteiro(const char *prompt)
{
    int v;
    while (true)
    {
        printf("%s", prompt);
        if (cin >> v)
            return v;
        cin.clear();
        cin.ignore(10000, '\n');
        printf(RED "  Entrada invalida. Digite um numero inteiro.\n" RST);
    }
}

double lerDouble(const char *prompt)
{
    double v;
    while (true)
    {
        printf("%s", prompt);
        if (cin >> v)
            return v;
        cin.clear();
        cin.ignore(10000, '\n');
        printf(RED "  Entrada invalida. Digite um numero.\n" RST);
    }
}

// Barra de progresso horizontal  [####....] pct%
void printBarra(double valor, double maximo, int largura, const char *cor)
{
    int cheio = (maximo > 0) ? (int)(valor / maximo * largura) : 0;
    if (cheio > largura)
        cheio = largura;
    printf("%s[", cor);
    for (int i = 0; i < largura; i++)
        putchar(i < cheio ? '#' : '.');
    printf("]" RST);
}

// Progresso de ingestão -reescreve a linha com \r
void printProgresso(unsigned int atual, unsigned int total)
{
    int pct = (total > 0) ? (int)((double)atual / total * 100) : 0;
    int blocos = pct * 40 / 100;
    printf("\r  " GRN "[");
    for (int i = 0; i < 40; i++)
        putchar(i < blocos ? '#' : '.');
    printf("]" RST " %3d%%  (%u amostras)", pct, atual);
    fflush(stdout);
}

// Cabeçalho do sistema
void printHeader()
{
    printLinha('=');
    printf(CYN BOLD
           "   SISTEMA DE ANALISE DE DETECCOES DE UVA\n"
           "   Estruturas de Dados aplicadas a Visao Computacional - 2026\n" RST);
    printLinha('=');
}

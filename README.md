# Calculadora Matricial - Teoria de Sistemas Lineares

Esta ferramenta é uma aplicação interativa desenvolvida em Python para auxiliar no estudo e na análise de sistemas lineares, com foco nas operações matriciais baseadas no Teorema de Cayley-Hamilton. A aplicação fornece um passo a passo detalhado (simbólico e numérico) para cada cálculo, promovendo um entendimento mais profundo dos conceitos matemáticos.

## Funcionalidades Principais

A aplicação possui uma interface gráfica intuitiva (Tkinter) dividida em três abas principais:

1.  **Potências de Matrizes (A^k):**
    *   Calcula $A^k$ para uma matriz $A$ de ordem $n$ e um expoente $k$.
    *   Utiliza a divisão polinomial simbólica de $\lambda^k$ pelo polinômio característico $p(\lambda)$.
    *   Exibe o polinômio característico, o resto $r(\lambda)$ e a expressão final $A^k = r(A)$.

2.  **Inversão de Matrizes (A^-1):**
    *   Calcula a inversa de uma matriz $A$ (se invertível).
    *   Verifica a singularidade da matriz (determinante zero).
    *   Exibe o polinômio característico e a fórmula da inversa derivada do Teorema de Cayley-Hamilton, com os coeficientes reais da matriz.
    *   Mostra o cálculo da soma interna e o resultado final da matriz inversa.

3.  **Resolução de EDOs Lineares (x' = Ax):**
    *   Resolve sistemas de equações diferenciais lineares homogêneas da forma $\dot{x} = Ax$.
    *   Calcula e exibe autovalores e autovetores.
    *   Fornece a matriz exponencial simbólica $e^{At}$ e a solução analítica $x(t) = e^{At}x_0$.
    *   Inclui um painel colapsável com gráficos interativos (Matplotlib) mostrando a evolução dos estados $x(t)$ ao longo do tempo.
    *   O painel de texto e o gráfico são dispostos lado a lado para melhor visualização.

## Pré-requisitos

*   Python 3.8 ou superior.

## Instalação

1.  **Clone o repositório** (ou baixe o arquivo `app_matrizes.py` diretamente):
    ```bash
    https://github.com/SafeMantella/TSL-matrix-calculator
    cd LinearSystemsTool-CayleyHamilton
    ```
    *(Nota: `https://github.com/SafeMantella/TSL-matrix-calculator` é um link ilustrativo. Substitua pelo link real do repositório se ele for disponibilizado no GitHub.)*

2.  **Instale as dependências** utilizando `pip` e o arquivo `requirements.txt`:
    ```bash
    pip install -r requirements.txt
    ```

## Como Executar

Após a instalação das dependências, execute o script principal:
```bash
python3 app_matrizes.py
```
A janela da aplicação gráfica será aberta.

## Uso

*   **Entrada da Matriz:** Em cada aba, primeiro insira a dimensão $n$ da matriz (entre 1 e 5). Clique em "Gerar Grid". Preencha os campos com os elementos da matriz. Você pode usar números inteiros, decimais, frações (ex: `1/2`) ou expressões simbólicas (ex: `sqrt(2)`).
*   **Cálculo:** Clique no botão "Calcular..." correspondente à funcionalidade desejada.
*   **Resultados:** O passo a passo e os resultados serão exibidos na área de texto à esquerda.
*   **Gráficos (EDOs):** Na aba de EDOs, um botão "Mostrar Gráficos ▶" permite expandir um painel com a visualização plotada.

## Contribuição

Contribuições são bem-vindas! Sinta-se à vontade para abrir "issues" para reportar bugs ou sugerir melhorias, e enviar "pull requests" com novas funcionalidades.

## Licença

Este projeto está licenciado sob a Licença MIT. Veja o arquivo `LICENSE` para mais detalhes.

## Autor

Pedro Arfux Pereira Cavalcante de Castro

## Agradecimentos

Agradeço à Universidade Federal de Mato Grosso do Sul (UFMS) pelo suporte institucional.

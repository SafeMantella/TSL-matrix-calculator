import tkinter as tk
from tkinter import ttk, messagebox, scrolledtext
import numpy as np
import sympy
from sympy import Matrix, eye, det, pprint, simplify, symbols
from sympy.abc import lamda as lambda_
from scipy.integrate import solve_ivp

# Matplotlib for plotting
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

class AppMatrizes:
    def __init__(self, root):
        self.root = root
        self.root.title("Calculadora Matricial - Teoria de Sistemas Lineares")
        self.root.geometry("900x850") # Aumentei um pouco a altura para caber o gráfico

        # Estilo
        style = ttk.Style()
        style.theme_use('clam')

        # Notebook (Abas)
        self.notebook = ttk.Notebook(root)
        self.notebook.pack(fill='both', expand=True, padx=10, pady=10)

        self.aba_potencias = ttk.Frame(self.notebook)
        self.aba_inversao = ttk.Frame(self.notebook)
        self.aba_edos = ttk.Frame(self.notebook)

        self.notebook.add(self.aba_potencias, text='Potências de Matrizes')
        self.notebook.add(self.aba_inversao, text='Inversão de Matrizes')
        self.notebook.add(self.aba_edos, text='Resolução de EDOs')

        # Configurar cada aba
        self.setup_potencias_tab()
        self.setup_inversao_tab()
        self.setup_edos_tab()

    def create_matrix_input_ui(self, parent, extra_setup_callback=None):
        """
        Cria a interface padrão de input de matriz (N e Grid).
        Retorna o frame do grid para referência.
        """
        frame_top = ttk.Frame(parent)
        frame_top.pack(fill='x', padx=10, pady=5)

        # Input N
        ttk.Label(frame_top, text="Tamanho da Matriz (n):").pack(side='left', padx=5)
        entry_n = ttk.Entry(frame_top, width=5)
        entry_n.pack(side='left', padx=5)
        entry_n.insert(0, "2")

        btn_gerar = ttk.Button(frame_top, text="Gerar Grid", 
                               command=lambda: self.gerar_grid(parent, entry_n, extra_setup_callback))
        btn_gerar.pack(side='left', padx=5)

        # Container para a matriz
        frame_grid = ttk.Frame(parent)
        frame_grid.pack(pady=10)

        # Armazenar referências no widget pai para acesso posterior
        parent.matrix_data = {
            'entry_n': entry_n,
            'frame_grid': frame_grid,
            'entries': []
        }
        
        return frame_top

    def create_output_area(self, parent):
        lbl = ttk.Label(parent, text="Resultados e Passos:")
        lbl.pack(anchor='w', padx=10, pady=(10, 0))
        txt = scrolledtext.ScrolledText(parent, height=15, width=80, font=('Courier', 10))
        txt.pack(fill='both', expand=True, padx=10, pady=5)
        return txt

    def gerar_grid(self, tab_frame, entry_n, extra_callback=None):
        try:
            n = int(entry_n.get())
            if n < 1 or n > 5:
                raise ValueError
        except:
            messagebox.showerror("Erro", "Por favor, insira um valor inteiro para n entre 1 e 5.")
            return

        data = tab_frame.matrix_data
        frame_grid = data['frame_grid']

        # Limpar grid anterior
        for widget in frame_grid.winfo_children():
            widget.destroy()

        data['entries'] = []

        # Criar novo grid
        for i in range(n):
            row_entries = []
            for j in range(n):
                e = ttk.Entry(frame_grid, width=8, justify='center')
                e.grid(row=i, column=j, padx=3, pady=3)
                row_entries.append(e)
            data['entries'].append(row_entries)

        # Callback para criar inputs extras dependentes de N (ex: vetor x0)
        if extra_callback:
            extra_callback(n)

    def ler_matriz(self, tab_frame):
        data = tab_frame.matrix_data
        entries = data['entries']
        if not entries:
            messagebox.showwarning("Aviso", "Gere a matriz primeiro.")
            return None

        n = len(entries)
        try:
            # Tenta converter inputs para SymPy (aceita inteiros, frações, sqrt, etc)
            matrix_list = []
            for i in range(n):
                row_vals = []
                for j in range(n):
                    val_str = entries[i][j].get()
                    if not val_str.strip():
                        raise ValueError(f"Célula ({i+1},{j+1}) vazia.")
                    row_vals.append(sympy.sympify(val_str))
                matrix_list.append(row_vals)
            return Matrix(matrix_list)
        except Exception as e:
            messagebox.showerror("Erro na Leitura", f"Erro ao ler a matriz: {e}")
            return None

    def print_to_output(self, text_widget, content, clear=False):
        if clear:
            text_widget.delete(1.0, tk.END)
        
        if isinstance(content, str):
            text_widget.insert(tk.END, content + "\n")
        else:
            # Use sympy.pretty directly to get string, avoiding pprint stream issues
            # Disable line wrapping with wrap_line=False and huge num_columns
            pretty_str = sympy.pretty(content, use_unicode=True, wrap_line=False, num_columns=10000)
            text_widget.insert(tk.END, pretty_str + "\n")
        
        text_widget.see(tk.END)

    # --- ABA 1: Potências ---
    def setup_potencias_tab(self):
        tab = self.aba_potencias
        
        self.create_matrix_input_ui(tab)
        
        # Input específico K
        frame_k = ttk.Frame(tab)
        frame_k.pack(pady=5)
        ttk.Label(frame_k, text="Expoente k:").pack(side='left')
        self.entry_k = ttk.Entry(frame_k, width=8)
        self.entry_k.pack(side='left', padx=5)
        
        btn_calc = ttk.Button(tab, text="Calcular Potência", command=self.calcular_potencia)
        btn_calc.pack(pady=5)
        
        tab.output = self.create_output_area(tab)

    def calcular_potencia(self):
        tab = self.aba_potencias
        M = self.ler_matriz(tab)
        if M is None: return

        try:
            k_val = int(self.entry_k.get())
        except:
            messagebox.showerror("Erro", "Insira um expoente k inteiro válido.")
            return

        out = tab.output
        self.print_to_output(out, f"--- Cálculo de A^{k_val} ---", clear=True)
        
        self.print_to_output(out, "Matriz A:")
        self.print_to_output(out, M)

        # Passo 1: Polinômio Característico
        lam = symbols('lambda')
        n = M.shape[0]
        p = M.charpoly(lam)
        p_expr = p.as_expr()
        
        self.print_to_output(out, f"\nPasso 1: Polinômio Característico p(λ) = det(A - λI):")
        self.print_to_output(out, p_expr)

        self.print_to_output(out, "\nPelo Teorema de Cayley-Hamilton, p(A) = 0.")
        
        # Passo 2: Divisão Polinomial
        self.print_to_output(out, f"\nPasso 2: Divisão de λ^{k_val} por p(λ):")
        self.print_to_output(out, f"Escrevemos λ^{k_val} = q(λ)*p(λ) + r(λ), onde o grau de r(λ) < n.")
        self.print_to_output(out, "Substituindo A na equação, temos A^k = q(A)*p(A) + r(A).")
        self.print_to_output(out, f"Como p(A) = 0, então A^{k_val} = r(A).")

        # Perform division
        term_k = lam**k_val
        quotient, remainder = sympy.div(term_k, p_expr, domain='QQ')
        
        self.print_to_output(out, f"\nResto da divisão r(λ):")
        self.print_to_output(out, remainder)
        
        # Passo 3: Calcular r(A)
        self.print_to_output(out, f"\nPasso 3: Calcular A^{k_val} substituindo A em r(λ):")
        
        # To calculate r(A), we need to extract coefficients of remainder
        r_poly = sympy.Poly(remainder, lam)
        r_coeffs = r_poly.all_coeffs()
        
        n = M.shape[0]
        res_matrix = sympy.zeros(n, n)
        
        calc_steps_str = []
        
        degree_r = r_poly.degree()
        
        # Check if remainder is effectively zero (empty list or zero val)
        if not r_coeffs or (len(r_coeffs) == 1 and r_coeffs[0] == 0):
             self.print_to_output(out, "O resto é 0, logo A^k = 0.")
             res_matrix = sympy.zeros(n,n)
        else:
            # r_coeffs are from highest degree down to constant
            # Example: 2*lam + 1 -> [2, 1] -> degree=1. 
            # i=0: coeff=2, power=1-0=1. i=1: coeff=1, power=1-1=0.
            
            for i, c_val in enumerate(r_coeffs):
                power = degree_r - i
                
                if power == 0:
                    term_mat = eye(n) * c_val
                    term_name = "I"
                else:
                    term_mat = (M**power) * c_val
                    term_name = f"A^{power}"
                
                res_matrix += term_mat
                
                # Formatting for display
                if c_val == 0: continue
                
                val_str = str(c_val)
                if i == 0:
                     calc_steps_str.append(f"{val_str}*{term_name}")
                else:
                    if c_val >= 0:
                        calc_steps_str.append(f"+ {val_str}*{term_name}")
                    else:
                        calc_steps_str.append(f"{val_str}*{term_name}")
            
            calc_str = " ".join(calc_steps_str)
            self.print_to_output(out, f"A^{k_val} = {calc_str}")
        
        self.print_to_output(out, "\nResultado Final:")
        self.print_to_output(out, res_matrix)

    # --- ABA 2: Inversão ---
    def setup_inversao_tab(self):
        tab = self.aba_inversao
        self.create_matrix_input_ui(tab)
        
        btn_calc = ttk.Button(tab, text="Calcular Inversa", command=self.calcular_inversa)
        btn_calc.pack(pady=10)
        
        tab.output = self.create_output_area(tab)

    def calcular_inversa(self):
        tab = self.aba_inversao
        M = self.ler_matriz(tab)
        if M is None: return

        out = tab.output
        self.print_to_output(out, "--- Cálculo da Inversa A^-1 ---", clear=True)
        self.print_to_output(out, "Matriz A:")
        self.print_to_output(out, M)

        det_val = M.det()
        
        if det_val == 0:
            self.print_to_output(out, f"\nDeterminante det(A) = {det_val}")
            self.print_to_output(out, "\nA matriz é singular (determinante zero), portanto não possui inversa.")
            return

        # Cayley-Hamilton Inversa
        lam = symbols('lambda')
        p = M.charpoly(lam)
        p_expr = p.as_expr()
        
        self.print_to_output(out, f"\nPasso 1: Polinômio Característico p(λ) = det(A - λI):")
        self.print_to_output(out, p_expr)

        coeffs = p.all_coeffs() # [1, c_{n-1}, ..., c_1, c_0]
        n = M.shape[0]
        c0 = coeffs[-1]

        self.print_to_output(out, "\nPasso 2: Aplicando o Teorema de Cayley-Hamilton:")
        self.print_to_output(out, f"Sabemos que p(A) = 0. Para o polinômio característico p(λ) = {p_expr}, temos:")
        
        # Build A^n + c_{n-1}A^{n-1} + ... + c_1 A + c_0 I = 0 dynamically
        cayley_hamilton_terms_list = []
        for i in range(n + 1):
            c_val = coeffs[i]
            power = n - i

            if c_val == 0: # Skip zero coefficients
                continue

            term_part = ""
            # Coefficient part
            if c_val == 1 and power > 0: # e.g., "A^2" instead of "1*A^2"
                coeff_str = ""
            elif c_val == -1 and power > 0: # e.g., "-A^2" instead of "-1*A^2"
                coeff_str = "-"
            else:
                coeff_str = str(c_val)

            # Matrix part
            if power == 0:
                matrix_str = "I"
            elif power == 1:
                matrix_str = "A"
            else:
                matrix_str = f"A^{power}"

            # Combine coeff and matrix part
            if coeff_str and coeff_str != "-": # e.g., "2*A^2"
                term_part = f"{coeff_str}*{matrix_str}"
            elif coeff_str == "-": # e.g., "-A^2"
                term_part = f"-{matrix_str}"
            else: # coeff_str is empty (meaning c_val is 1)
                term_part = matrix_str
            
            cayley_hamilton_terms_list.append(term_part)
        
        # Join the terms, handling signs for readability
        final_ch_eq_parts = []
        for term in cayley_hamilton_terms_list:
            if term.startswith('-') or not final_ch_eq_parts: # If negative or the first term
                final_ch_eq_parts.append(term)
            else: # If positive and not the first term
                final_ch_eq_parts.append(f"+ {term}")
        
        cayley_hamilton_equation_str = " ".join(final_ch_eq_parts) + " = 0"
        self.print_to_output(out, cayley_hamilton_equation_str)        
        
        self.print_to_output(out, "\nMultiplicando por A^-1 e isolando o termo, obtemos a fórmula da inversa:")
        
        # Construir string da fórmula com valores
        terms = []
        for i in range(n):
            c_val = coeffs[i]
            power = n - 1 - i
            
            term_str = ""
            
            # Sinal
            if i > 0 and c_val >= 0:
                term_str += "+ "
            elif c_val < 0:
                term_str += "- " 
            
            val_abs = abs(c_val)
            
            if power == 0:
                base = "I"
            else:
                base = f"A^{power}"
            
            if val_abs == 1:
                term_str += f"{base}"
            else:
                term_str += f"{val_abs}*{base}"
            
            terms.append(term_str)
            
        # Juntar
        # Substituir o sinal de - separado se necessário
        terms_simple = []
        for i in range(n):
            c_val = coeffs[i]
            power = n - 1 - i
            
            if power == 0:
                base = "I"
            else:
                base = f"A^{power}"
                
            if i == 0:
                # Primeiro termo (geralmente 1*A^n-1)
                terms_simple.append(f"{c_val}*{base}")
            else:
                if c_val >= 0:
                    terms_simple.append(f"+ {c_val}*{base}")
                else:
                    terms_simple.append(f"{c_val}*{base}") 
        
        poly_str = " ".join(terms_simple)
        
        self.print_to_output(out, f"A^-1 = (-1/{c0}) * ({poly_str})")
        
        self.print_to_output(out, f"\nCoeficientes identificados:")
        for i, c in enumerate(coeffs):
            power = n - i
            self.print_to_output(out, f"Coeficiente de λ^{power} (c_{power}): {c}")
        
        self.print_to_output(out, f"\nTermo constante c_0 = {c0} (Note que c_0 = (-1)^n * det(A))")
        
        # Calcular a soma interna: (A^{n-1} + c_{n-1}A^{n-2} + ... + c_1 I)
        self.print_to_output(out, "\nPasso 3: Calculando a soma interna S = (A^{n-1} + ... + c_1 I):")
        
        S = sympy.zeros(n, n)
        
        formula_str_parts = []
        
        for i in range(n): # i vai de 0 a n-1. coeffs[i] corresponde a c_{n-i}
            c_val = coeffs[i]
            power_of_A = n - 1 - i
            
            term_matrix = (M ** power_of_A) * c_val
            S += term_matrix
            
            term_str = f"({c_val}) * A^{power_of_A}"
            if power_of_A == 0:
                term_str = f"({c_val}) * I"
            formula_str_parts.append(term_str)

        formula_str = " + ".join(formula_str_parts)
        self.print_to_output(out, f"S = {formula_str}")
        self.print_to_output(out, "Matriz S resultante:")
        self.print_to_output(out, S)

        self.print_to_output(out, f"\nPasso 4: Multiplicando por -1/c_0 = -1/({c0}):")
        
        inv = S * (-1/c0)
        self.print_to_output(out, "Resultado Final A^-1:")
        self.print_to_output(out, inv)

    # --- ABA 3: EDOs ---
    def setup_edos_tab(self):
        tab = self.aba_edos
        
        # Callback para gerar inputs do vetor x0 quando o grid for gerado
        self.create_matrix_input_ui(tab, extra_setup_callback=self.setup_edos_x0)
        
        # Frame para x0 (será preenchido pelo callback)
        self.frame_x0_container = ttk.LabelFrame(tab, text="Vetor Inicial x(0)")
        self.frame_x0_container.pack(pady=5, fill='x', padx=10)
        
        # Frame para tempo
        frame_t = ttk.Frame(tab)
        frame_t.pack(pady=5)
        ttk.Label(frame_t, text="Tempo final t:").pack(side='left')
        self.entry_t = ttk.Entry(frame_t, width=8)
        self.entry_t.insert(0, "10")
        self.entry_t.pack(side='left', padx=5)

        btn_solve = ttk.Button(tab, text="Resolver EDO", command=self.resolver_edo)
        btn_solve.pack(pady=5)

        # Botão Toggle Plot
        self.btn_toggle_plot = ttk.Button(tab, text="Mostrar Gráficos ▶", command=self.toggle_plot)
        self.btn_toggle_plot.pack(pady=5)
        
        # Container Principal para Resultados (Lado a Lado)
        self.results_container = ttk.Frame(tab)
        self.results_container.pack(fill='both', expand=True, padx=5, pady=5)

        # --- Interface de Texto (Esquerda) ---
        # Precisamos recriar a area de output sem o pack padrão do helper se quisermos controlar aqui,
        # ou passar o container. O helper create_output_area faz pack fill both.
        # Vamos usar um frame wrapper para a esquerda.
        self.frame_text_wrapper = ttk.Frame(self.results_container)
        self.frame_text_wrapper.pack(side='left', fill='both', expand=True)
        
        # Cria o output dentro do wrapper
        tab.output = self.create_output_area(self.frame_text_wrapper)

        # --- Interface de Gráficos (Direita) ---
        self.frame_plot = ttk.Frame(self.results_container)
        # O pack será feito no toggle
        
        # Configuração do Matplotlib
        self.fig, self.ax = plt.subplots(figsize=(5, 4), dpi=100) # Ajustei figsize
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.frame_plot)
        self.canvas_widget = self.canvas.get_tk_widget()
        self.canvas_widget.pack(fill='both', expand=True)

    def toggle_plot(self):
        if self.frame_plot.winfo_viewable():
            self.frame_plot.pack_forget()
            self.btn_toggle_plot.config(text="Mostrar Gráficos ▶")
        else:
            self.frame_plot.pack(side='right', fill='both', expand=True, padx=(5, 0))
            self.btn_toggle_plot.config(text="Ocultar Gráficos ◀")

    def setup_edos_x0(self, n):
        # Limpar container x0
        for w in self.frame_x0_container.winfo_children():
            w.destroy()
            
        self.aba_edos.x0_entries = []
        for i in range(n):
            e = ttk.Entry(self.frame_x0_container, width=6)
            e.pack(side='left', padx=2, pady=5)
            e.insert(0, "0") # default
            self.aba_edos.x0_entries.append(e)

    def resolver_edo(self):
        tab = self.aba_edos
        M = self.ler_matriz(tab)
        if M is None: return

        out = tab.output
        self.print_to_output(out, "--- Resolução de x' = Ax ---", clear=True)
        
        # Ler x0
        try:
            x0 = [float(e.get()) for e in tab.x0_entries]
        except:
            messagebox.showerror("Erro", "Valores inválidos em x(0). Use números reais.")
            return

        try:
            tf = float(self.entry_t.get())
        except:
            messagebox.showerror("Erro", "Tempo inválido.")
            return

        self.print_to_output(out, "Sistema:")
        self.print_to_output(out, M)
        self.print_to_output(out, f"\nVetor inicial x(0): {x0}")
        self.print_to_output(out, f"Intervalo de tempo: [0, {tf}]")

        # Passo 1: Autovalores/Autovetores
        self.print_to_output(out, "\nPasso 1: Autovalores e Autovetores:")
        try:
            eigen_info = M.eigenvects()
            for val, mult, vecs in eigen_info:
                self.print_to_output(out, f"Autovalor: {val} (multiplicidade {mult})")
                for v in vecs:
                     self.print_to_output(out, f"  Autovetor: {list(v)}")
        except Exception as e:
            self.print_to_output(out, f"Erro ao calcular autovalores: {e}")

        # Passo 2: Matriz Fundamental / Exponencial
        self.print_to_output(out, "\nPasso 2: Matriz de Transição de Estados Φ(t) = e^{At}:")
        t_sym = symbols('t')
        expAt = None
        try:
            # Try compute exp(At) explicitly
            expAt = (M * t_sym).exp()
            self.print_to_output(out, expAt)
        except:
             self.print_to_output(out, "Cálculo simbólico complexo demais para exibir em grade ou erro no SymPy.")
        
        # Passo 3: Solução Geral
        self.print_to_output(out, "\nPasso 3: Solução x(t) = e^{At} * x(0):")
        if expAt is not None:
            try:
                x0_mat = Matrix(x0)
                sol_sym = expAt * x0_mat
                self.print_to_output(out, sol_sym)
            except Exception as e:
                 self.print_to_output(out, f"Erro na multiplicação simbólica: {e}")

        # Solução Numérica (SciPy)
        self.print_to_output(out, "\n--- Solução Numérica Final (t = final) ---")
        
        # Converter matriz sympy para numpy float
        try:
            A_np = np.array(M.tolist(), dtype=float)
        except:
            self.print_to_output(out, "Aviso: Matriz contém símbolos não numéricos, impossível resolver numericamente.")
            return

        def sistema(t, x):
            return A_np @ x

        sol = solve_ivp(sistema, [0, tf], x0, t_eval=np.linspace(0, tf, 200))
        
        # Mostrar resultado final
        xf = sol.y[:, -1]
        self.print_to_output(out, f"x({tf}) = \n{xf}")

        # --- Plotagem dos Gráficos ---
        try:
            self.ax.clear()
            # sol.y shape é (n, pontos_tempo)
            for i in range(sol.y.shape[0]):
                self.ax.plot(sol.t, sol.y[i], label=f'$x_{i+1}(t)$', linewidth=2)
            
            self.ax.set_title("Evolução dos Estados no Tempo")
            self.ax.set_xlabel("Tempo (t)")
            self.ax.set_ylabel("Estados x(t)")
            self.ax.grid(True, linestyle='--', alpha=0.7)
            self.ax.legend()
            
            self.canvas.draw()
            
            # Abre automaticamente o painel de gráficos se estiver fechado
            if not self.frame_plot.winfo_viewable():
                self.toggle_plot()
                
        except Exception as e:
            self.print_to_output(out, f"\nErro ao gerar gráfico: {e}")

if __name__ == "__main__":
    root = tk.Tk()
    app = AppMatrizes(root)
    root.mainloop()
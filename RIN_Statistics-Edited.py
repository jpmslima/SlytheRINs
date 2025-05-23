#!/usr/bin/env python
# coding: utf-8

# # 🧬 Análise de Redes com Estatísticas e Visualizações Completas

# In[1]:


get_ipython().system('pip install -U kaleido')


# In[4]:


# Importação das bibliotecas necessárias
import os
import glob
import pandas as pd
import numpy as np
import networkx as nx
import plotly.graph_objects as go
import plotly.express as px
from pandas.errors import EmptyDataError
import scipy.stats as st

print("SlytheRINs - Network Analysis Application")
print("Analyzing Residue Interaction Networks (RINs)")
print("---------------------------------------------------------")

# ======================================================
# Parte 1: Carregando os arquivos .edges.txt da pasta atual
# ======================================================
# Busca por arquivos que terminem com ".edges.txt" na pasta corrente
edge_file_list = glob.glob("*.edges.txt")

if not edge_file_list:
    print("Nenhum arquivo '.edges.txt' encontrado na pasta atual.")
    edge_dataframes = []
else:
    edge_dataframes = []
    for edge_file in edge_file_list:
        try:
            df = pd.read_csv(edge_file, sep='\t')
            edge_dataframes.append(df)
            print(f"Arquivo carregado: {edge_file}")
        except EmptyDataError:
            print(f"O arquivo {edge_file} está vazio ou não possui colunas a serem interpretadas.")
        except Exception as e:
            print(f"Ocorreu um erro ao ler o arquivo {edge_file}: {e}")

# ======================================================
# Parte 2: Processamento dos dados de rede e criação dos gráficos
# ======================================================
print("\n*** First Part: Network Analysis and Comparison Results ***")
if edge_dataframes:
    # Inicializa listas para armazenar parâmetros de cada rede
    all_degrees = []
    all_betweenness = []
    all_clustering = []
    all_triangles = []
    assortativity_values = []

    for df in edge_dataframes:
        # Cria grafo usando as colunas 'NodeId1' e 'NodeId2'
        G = nx.from_pandas_edgelist(df, 'NodeId1', 'NodeId2')
        degree = dict(G.degree())
        betweenness = nx.betweenness_centrality(G)
        clustering = nx.clustering(G)
        triangles = nx.triangles(G)
        assortativity = nx.degree_assortativity_coefficient(G)

        all_degrees.append(degree)
        all_betweenness.append(betweenness)
        all_clustering.append(clustering)
        all_triangles.append(triangles)
        assortativity_values.append(assortativity)

    # Cálculo da variância dos valores de grau, betweenness e clustering
    all_degree_values = pd.DataFrame(all_degrees).fillna(0)
    degree_variance = all_degree_values.var(axis=0)

    all_betweenness_values = pd.DataFrame(all_betweenness).fillna(0)
    betweenness_variance = all_betweenness_values.var(axis=0)

    all_clustering_values = pd.DataFrame(all_clustering).fillna(0)
    clustering_variance = all_clustering_values.var(axis=0)

    # Cálculo da média e desvio padrão dos valores de grau
    mean_degree = all_degree_values.mean(axis=0)
    std_degree = all_degree_values.std(axis=0)
    neg_std_degree = -std_degree

    # Ordenação dos resíduos – supõe-se que os índices (NodeId1) estejam no formato "Chain:Position:res:AA"
    residue_numbers = degree_variance.index.tolist()
    sorted_residue_numbers = sorted(residue_numbers, key=lambda x: int(x.split(':')[1]))
    formatted_labels = [f'{x.split(":")[1]}:{x.split(":")[3]}' for x in sorted_residue_numbers]

    # Processamento dos dados dos triângulos
    triangle_counts_df = pd.DataFrame(all_triangles).fillna(0)
    mean_triangles = triangle_counts_df.mean(axis=0)
    std_triangles = triangle_counts_df.std(axis=0)

    # ---------------------------------------------------------------------
    # Salvando os dados utilizados para os gráficos em arquivos TSV
    # ---------------------------------------------------------------------
    # 1. Dados de variância do grau
    df_degree_variance = pd.DataFrame({
        'Residue': sorted_residue_numbers,
        'Formatted_Label': formatted_labels,
        'Degree_Variance': degree_variance.loc[sorted_residue_numbers].values
    })
    df_degree_variance.to_csv("degree_variance.tsv", sep="\t", index=False)
    print("Dados de variância do grau salvos em 'degree_variance.tsv'.")

    # 2. Dados de variância do betweenness
    df_betweenness_variance = pd.DataFrame({
        'Residue': sorted_residue_numbers,
        'Formatted_Label': formatted_labels,
        'Betweenness_Variance': betweenness_variance.loc[sorted_residue_numbers].values
    })
    df_betweenness_variance.to_csv("betweenness_variance.tsv", sep="\t", index=False)
    print("Dados de variância do betweenness salvos em 'betweenness_variance.tsv'.")

    # 3. Dados de variância do clustering
    df_clustering_variance = pd.DataFrame({
        'Residue': sorted_residue_numbers,
        'Formatted_Label': formatted_labels,
        'Clustering_Variance': clustering_variance.loc[sorted_residue_numbers].values
    })
    df_clustering_variance.to_csv("clustering_variance.tsv", sep="\t", index=False)
    print("Dados de variância do clustering salvos em 'clustering_variance.tsv'.")

    # 4. Dados de média e STD do grau
    df_mean_std_degree = pd.DataFrame({
        'Residue': sorted_residue_numbers,
        'Formatted_Label': formatted_labels,
        'Mean_Degree': mean_degree.loc[sorted_residue_numbers].values,
        'STD_Degree': std_degree.loc[sorted_residue_numbers].values
    })
    df_mean_std_degree.to_csv("mean_std_degree.tsv", sep="\t", index=False)
    print("Dados de média e desvio padrão do grau salvos em 'mean_std_degree.tsv'.")

    # 5. Dados dos triângulos (média e STD)
    df_triangle_counts = pd.DataFrame({
        'Residue': sorted_residue_numbers,
        'Formatted_Label': formatted_labels,
        'Mean_Triangles': mean_triangles.loc[sorted_residue_numbers].values,
        'STD_Triangles': std_triangles.loc[sorted_residue_numbers].values
    })
    df_triangle_counts.to_csv("triangle_counts.tsv", sep="\t", index=False)
    print("Dados dos triângulos salvos em 'triangle_counts.tsv'.")

    # 6. Dados de assortatividade (um valor por grafo)
    df_assortativity = pd.DataFrame({
        'Graph_Index': list(range(1, len(assortativity_values) + 1)),
        'Assortativity': assortativity_values
    })
    df_assortativity.to_csv("assortativity.tsv", sep="\t", index=False)
    print("Dados de assortatividade salvos em 'assortativity.tsv'.")
    # ---------------------------------------------------------------------

    # ---------------------------------------------------------------------
    # Criação dos gráficos com Plotly e salvando em HTML e PNG
    # ---------------------------------------------------------------------
    # Função auxiliar para salvar os gráficos em HTML e PNG
    def save_graph(fig, base_filename):
        fig.write_html(f"{base_filename}.html")
        print(f"Gráfico salvo em '{base_filename}.html'.")
        # Salva como PNG: certifique-se de ter kaleido instalado.
        fig.write_image(f"{base_filename}.png")
        print(f"Gráfico salvo em '{base_filename}.png'.")

    # Gráfico: Variância do grau
    fig_degree = go.Figure()
    fig_degree.add_trace(go.Scatter(
        x=sorted_residue_numbers,
        y=degree_variance.loc[sorted_residue_numbers].values,
        mode='lines+markers'
    ))
    fig_degree.update_layout(
        title='Variance of Degree Values Across Residues',
        xaxis_title='Residue Position and Amino Acid',
        yaxis_title='Variance of Degree Values',
        xaxis_tickangle=-90
    )
    fig_degree.show()
    save_graph(fig_degree, "graph_degree")

    # Gráfico: Variância do betweenness
    fig_betweenness = go.Figure()
    fig_betweenness.add_trace(go.Scatter(
        x=sorted_residue_numbers,
        y=betweenness_variance.loc[sorted_residue_numbers].values,
        mode='lines+markers'
    ))
    fig_betweenness.update_layout(
        title='Variance of Betweenness Centrality Across Residues',
        xaxis_title='Residue Position and Amino Acid',
        yaxis_title='Variance of Betweenness Centrality',
        xaxis_tickangle=-90
    )
    fig_betweenness.show()
    save_graph(fig_betweenness, "graph_betweenness")

    # Gráfico: Variância do clustering
    fig_clustering = go.Figure()
    fig_clustering.add_trace(go.Scatter(
        x=sorted_residue_numbers,
        y=clustering_variance.loc[sorted_residue_numbers].values,
        mode='lines+markers'
    ))
    fig_clustering.update_layout(
        title='Variance of Clustering Coefficient Across Residues',
        xaxis_title='Residue Position and Amino Acid',
        yaxis_title='Variance of Clustering Coefficient',
        xaxis_tickangle=-90
    )
    fig_clustering.show()
    save_graph(fig_clustering, "graph_clustering")

    # Gráfico: Triângulos (média e STD)
    fig_triangles = go.Figure()
    fig_triangles.add_trace(go.Scatter(
        x=sorted_residue_numbers,
        y=mean_triangles.loc[sorted_residue_numbers].values,
        mode='lines+markers',
        name='Mean Triangles'
    ))
    fig_triangles.add_trace(go.Scatter(
        x=sorted_residue_numbers,
        y=(mean_triangles - std_triangles).loc[sorted_residue_numbers].values,
        mode='lines',
        fill=None,
        line_color='lightblue',
        name='Lower Bound'
    ))
    fig_triangles.add_trace(go.Scatter(
        x=sorted_residue_numbers,
        y=(mean_triangles + std_triangles).loc[sorted_residue_numbers].values,
        mode='lines',
        fill='tonexty',
        line_color='lightblue',
        name='Upper Bound'
    ))
    fig_triangles.update_layout(
        title='Mean and Standard Deviation of Triangle Counts per Residue',
        xaxis_title='Residue',
        yaxis_title='Number of Triangles'
    )
    fig_triangles.show()
    save_graph(fig_triangles, "graph_triangles")

    # Gráfico: Média e STD do grau
    fig_degree_sd = go.Figure()
    fig_degree_sd.add_trace(go.Scatter(
        x=sorted_residue_numbers,
        y=mean_degree.loc[sorted_residue_numbers].values,
        mode='lines+markers',
        name='Mean Degree'
    ))
    fig_degree_sd.add_trace(go.Scatter(
        x=sorted_residue_numbers,
        y=(mean_degree - std_degree).loc[sorted_residue_numbers].values,
        mode='lines',
        fill=None,
        line_color='lightblue',
        name='Lower Bound'
    ))
    fig_degree_sd.add_trace(go.Scatter(
        x=sorted_residue_numbers,
        y=(mean_degree + std_degree).loc[sorted_residue_numbers].values,
        mode='lines',
        fill='tonexty',
        line_color='lightblue',
        name='Upper Bound'
    ))
    fig_degree_sd.update_layout(
        title='Mean Degree Value per Residue Position',
        xaxis_title='Residue Position and Amino Acid',
        yaxis_title='Degree',
        xaxis_tickangle=-90
    )
    fig_degree_sd.show()
    save_graph(fig_degree_sd, "graph_mean_std_degree")

    # Gráfico: Assortatividade dos grafos
    fig_assortativity = go.Figure()
    fig_assortativity.add_trace(go.Scatter(
        x=list(range(1, len(assortativity_values) + 1)),
        y=assortativity_values,
        mode='lines+markers',
        name='Assortativity Coefficient'
    ))
    fig_assortativity.update_layout(
        title='Assortativity Coefficient for Each Graph',
        xaxis_title='Graph Index',
        yaxis_title='Assortativity Coefficient'
    )
    fig_assortativity.show()
    save_graph(fig_assortativity, "graph_assortativity")

    # Gráfico: STD e -STD do grau com regiões sombreadas
    fig_dsd = go.Figure()
    fig_dsd.add_trace(go.Scatter(
        x=sorted_residue_numbers,
        y=std_degree.loc[sorted_residue_numbers].values,
        mode='lines+markers',
        name='Standard Deviation'
    ))
    fig_dsd.add_trace(go.Scatter(
        x=sorted_residue_numbers,
        y=neg_std_degree.loc[sorted_residue_numbers].values,
        mode='lines+markers',
        name='Negative Standard Deviation'
    ))
    fig_dsd.add_shape(type='rect', x0=sorted_residue_numbers[0], x1=sorted_residue_numbers[-1],
                      y0=-0.5, y1=0.5, fillcolor='lightgrey', opacity=0.5, layer='below', line_width=0)
    fig_dsd.add_shape(type='rect', x0=sorted_residue_numbers[0], x1=sorted_residue_numbers[-1],
                      y0=-1.0, y1=1.0, fillcolor='lightgrey', opacity=0.3, layer='below', line_width=0)
    fig_dsd.update_layout(
        title='Standard Deviation and Negative Standard Deviation of Degree Values',
        xaxis_title='Residue Position',
        yaxis_title='Relative Standard Deviation of Degree'
    )
    fig_dsd.show()
    save_graph(fig_dsd, "graph_std_negstd_degree")

else:
    print("Nenhum arquivo '.edges.txt' foi carregado. Verifique a pasta atual.")

# ---------------------------------------------------------------------
# Adicional: Teste T pareado para os valores de degree por resíduo
# ---------------------------------------------------------------------
paired_t_p_values = []

for residue in sorted_residue_numbers:
    # Valores do grau para o resíduo "residue" em cada rede
    x = all_degree_values[residue].values

    # Para cada rede, calcule a média dos graus EXCLUINDO o resíduo em questão
    overall_means_excluding = []
    for idx, row in all_degree_values.iterrows():
        total = row.sum() - row[residue]
        count = len(row) - 1
        overall_mean_excluding = total / count if count > 0 else 0
        overall_means_excluding.append(overall_mean_excluding)
    overall_means_excluding = pd.Series(overall_means_excluding, index=all_degree_values.index)

    # Realiza o teste t pareado: H0 => diferença média = 0 entre os graus do resíduo e a média das demais posições
    t_stat, p_val_t = st.ttest_rel(x, overall_means_excluding)
    paired_t_p_values.append(p_val_t)

# Cria um DataFrame com os resultados: resíduo, label formatado, variação (já calculada) e p-value do teste t pareado
df_significance = pd.DataFrame({
    "Residue": sorted_residue_numbers,
    "Formatted_Label": formatted_labels,
    "Degree_Variance": degree_variance.loc[sorted_residue_numbers].values,
    "Paired_T_p_value": paired_t_p_values,
})

df_significance.to_csv("residue_significance.tsv", sep="\t", index=False)
print("Arquivo 'residue_significance.tsv' com variação e p-values (teste t pareado) gerado.")

# Geração do gráfico dos p-values para cada resíduo (apenas o teste t pareado)
fig_pvalues = go.Figure()
fig_pvalues.add_trace(go.Scatter(
    x=sorted_residue_numbers,
    y=paired_t_p_values,
    mode='lines+markers',
    name="Paired T-Test P-Value"
))

# Linha horizontal para o nível de significância 0.05 (95%)
fig_pvalues.add_shape(
    dict(
        type="line",
        x0=sorted_residue_numbers[0],
        x1=sorted_residue_numbers[-1],
        y0=0.05,
        y1=0.05,
        line=dict(color="red", dash="dash")
    )
)
fig_pvalues.update_layout(
    title="P-Values por Resíduo (Teste T Pareado)",
    xaxis_title="Resíduo (ordenado)",
    yaxis_title="P-Value",
    xaxis_tickangle=-90
)
fig_pvalues.show()
save_graph(fig_pvalues, "p_values_graph")



# ---------------------------------------------------------------------
# Adicional: Teste T pareado para assortatividade
# ---------------------------------------------------------------------
# (Re)definição de 'assort_df' caso a célula original não tenha sido executada
assort_df = pd.Series(assortatividade, name='Assortatividade')

from scipy.stats import ttest_1samp
import plotly.graph_objects as go

# --- Teste t de uma amostra para assortatividade contra popmean=0 ---
stat_assort, p_assort = ttest_1samp(assort_df, popmean=0)
print(f"Assortatividade média = {assort_df.mean():.4f}  |  t = {stat_assort:.3f}  |  p = {p_assort:.3e}")

# --- Gráfico da assortatividade por modelo, com anotação de p<0.05 ---
signif = p_assort < 0.05
fig_assort = go.Figure()
fig_assort.add_trace(go.Scatter(
    x=list(range(len(assort_df))),
    y=assort_df,
    mode='lines+markers',
    name='Assortatividade'
))
# linha horizontal em zero
fig_assort.add_hline(y=0, line_dash="dash", annotation_text="zero", annotation_position="bottom right")

# Anotação de significância e p-valor
texto = f"{'* ' if signif else ''}p = {p_assort:.3e}"
fig_assort.add_annotation(
    xref="paper", yref="paper",
    x=0.98, y=0.98,
    text=texto,
    showarrow=False,
    font=dict(size=12)
)

fig_assort.update_layout(
    title="📊 Assortatividade por Modelo" + ("  *p<0.05*" if signif else ""),
    xaxis_title="Modelo",
    yaxis_title="Assortatividade",
    margin=dict(t=80)
)

# Salvar e exibir
fig_assort.write_html("Assortatividade_teste_t.html")
fig_assort.write_image("Assortatividade_teste_t.png")
print("Gráfico salvo como Assortatividade_teste_t.html e .png")
fig_assort.show()

import scipy.stats as st
import numpy as np
import pandas as pd
import plotly.graph_objects as go

# --- Make sure 'all_degree_values' DataFrame is available ---
# Example: all_degree_values = pd.DataFrame(all_degrees).fillna(0)

print("\n*** Additional Part: Testing Significance of Degree Variation per Residue ***")

variation_p_values = []
residue_list = all_degree_values.columns.tolist()
sorted_residue_numbers = sorted(residue_list, key=lambda x: int(x.split(':')[1]))
formatted_labels = [f'{x.split(":")[1]}:{x.split(":")[3]}' for x in sorted_residue_numbers]


N = len(all_degree_values)  # Number of conformations/networks

if N <= 1:
    print("Need more than one network to calculate variation significance.")
else:
    for residue in sorted_residue_numbers:
        degrees = all_degree_values[residue].values
        mean_deg = np.mean(degrees)
        var_deg = np.var(degrees, ddof=1)  # Use sample variance (ddof=1)

        # Handle cases with no variation or zero mean
        if var_deg == 0 or mean_deg == 0:
            p_value = 1.0  # No variation or zero mean -> not significantly variable
        else:
            # Chi-squared statistic for VMR test
            chi2_stat = (N - 1) * var_deg / mean_deg

            # P-value from Chi-squared distribution with N-1 degrees of freedom
            # We calculate a two-tailed p-value to test if variance is significantly
            # *different* (either larger or smaller) from the mean.
            cdf_val = st.chi2.cdf(chi2_stat, N - 1)
            p_value = 2 * min(cdf_val, 1 - cdf_val)

        variation_p_values.append(p_value)

    # Create a DataFrame with the results
    df_variation_pvals = pd.DataFrame({
        "Residue": sorted_residue_numbers,
        "Formatted_Label": formatted_labels,
        "Degree_Variation_P_Value": variation_p_values,
    })

    df_variation_pvals.to_csv("degree_variation_significance.tsv", sep="\t", index=False)
    print("Arquivo 'degree_variation_significance.tsv' com p-values para variação do grau gerado.")

    # Generate the plot
    fig_var_pvals = go.Figure()
    fig_var_pvals.add_trace(go.Scatter(
        x=sorted_residue_numbers,
        y=variation_p_values,
        mode='lines+markers',
        name="Degree Variation P-Value (vs. Poisson)"
    ))

    # Add 0.05 significance line
    fig_var_pvals.add_shape(
        dict(
            type="line",
            x0=sorted_residue_numbers[0],
            x1=sorted_residue_numbers[-1],
            y0=0.05,
            y1=0.05,
            line=dict(color="red", dash="dash")
        )
    )
    fig_var_pvals.update_layout(
        title="P-Values for Degree Variation per Residue (vs. Poisson)",
        xaxis_title="Residue (ordenado)",
        yaxis_title="P-Value",
        xaxis_tickangle=-90
    )
    fig_var_pvals.show()
    # Assuming 'save_graph' function is defined as in your script
    save_graph(fig_var_pvals, "p_values_variation_graph")

# ======================================================
# Parte 3: Análise das Interações Químicas
# ======================================================
print("\n*** Second Part: Analysis of the Chemical Interactions ***")
if edge_dataframes:
    # Conta as interações por tipo e por resíduo
    interaction_counts = {}
    for df in edge_dataframes:
        if 'Interaction' in df.columns and 'NodeId1' in df.columns:
            counts = df.groupby(['Interaction', 'NodeId1']).size().reset_index(name='Count')
            for interaction in counts['Interaction'].unique():
                if interaction not in interaction_counts:
                    interaction_counts[interaction] = []
                interaction_counts[interaction].append(counts[counts['Interaction'] == interaction])
        else:
            print("Os dados não contêm as colunas requeridas 'Interaction' e 'NodeId1'.")

    # Calcula a média e o STD para cada tipo de interação por resíduo
    interaction_stats_counts = {}
    for interaction, data_list in interaction_counts.items():
        combined_data = pd.concat(data_list)
        stats = combined_data.groupby('NodeId1')['Count'].agg(['mean', 'std']).reset_index()
        # Formatação de NodeId1
        def format_node_id(node_id):
            parts = node_id.split(':')
            return f'{parts[1]}:{parts[3]}'
        stats['FormattedNodeId1'] = stats['NodeId1'].apply(format_node_id)
        interaction_stats_counts[interaction] = stats
        # Salva os dados de cada interação em arquivo TSV
        file_name = f"interaction_{interaction}_counts.tsv"
        stats.to_csv(file_name, sep="\t", index=False)
        print(f"Dados de interações de '{interaction}' salvos em '{file_name}'.")

    # Função para criar gráfico de contagem de interações
    def create_interaction_count_plot(category, stats):
        figInt = go.Figure()
        figInt.add_trace(go.Scatter(
            x=stats['FormattedNodeId1'],
            y=stats['mean'],
            error_y=dict(type='data', array=stats['std']),
            mode='markers',
            marker=dict(size=10),
            name=category
        ))
        figInt.update_layout(
            title=f'{category} Interaction Counts',
            xaxis_title='Residue Position',
            yaxis_title='Mean Count',
            hovermode='closest'
        )
        return figInt

    # Gerando os gráficos para cada tipo de interação
    for category, stats in interaction_stats_counts.items():
        print(f"Plotando contagem de interações para {category}")
        figInt = create_interaction_count_plot(category, stats)
        figInt.show()
        base_filename = f"graph_interaction_{category}"
        save_graph(figInt, base_filename)
else:
    print("Não há dados de edges para análise das interações.")

# ======================================================
# Parte 4: Integração com AlphaMissense
# ======================================================
print("\n*** Third Part: Integration with AlphaMissense predicted mutation effects ***")
# Para a integração com AlphaMissense, o código busca arquivos TSV que contenham "AlphaMissense" no nome.
am_file_list = glob.glob("*.tsv")
am_file_list = [f for f in am_file_list if "AlphaMissense" in f]

if not am_file_list:
    print("Nenhum arquivo AlphaMissense encontrado na pasta atual.")
    am_df = pd.DataFrame()
else:
    am_df = pd.DataFrame()
    # Se houver mais de um arquivo, concatena-os
    for am_file in am_file_list:
        try:
            temp_df = pd.read_csv(am_file, sep='\t')
            am_df = pd.concat([am_df, temp_df], ignore_index=True)
            print(f"Arquivo AlphaMissense carregado: {am_file}")
        except EmptyDataError:
            print(f"O arquivo {am_file} está vazio ou não possui colunas para leitura.")
        except Exception as e:
            print(f"Ocorreu um erro ao ler o arquivo {am_file}: {e}")

if not am_df.empty:
    if am_df.shape[1] > 1:
        # Supondo que a coluna de posição seja 'Mutation_Pos'; se não, usa a segunda coluna
        if 'Mutation_Pos' in am_df.columns:
            am_df['Position'] = am_df['Mutation_Pos'].astype(int)
        else:
            am_df['Position'] = am_df.iloc[:, 1].astype(int)

        # Salva os dados completos do AlphaMissense
        am_df.to_csv("AlphaMissense_data.tsv", sep="\t", index=False)
        print("Dados do AlphaMissense salvos em 'AlphaMissense_data.tsv'.")

        # Conta as ocorrências de cada classificação por posição
        if 'Classification' in am_df.columns:
            classification_counts = am_df.groupby(['Position', 'Classification']).size().unstack(fill_value=0)
            classification_counts.to_csv("AlphaMissense_classification_counts.tsv", sep="\t")
            print("Contagens de classificação salvos em 'AlphaMissense_classification_counts.tsv'.")
        else:
            print("A coluna 'Classification' não foi encontrada no arquivo AlphaMissense.")

        # Define o score a partir de uma coluna (por ex.: 'Score' ou similar)
        if 'Score' in am_df.columns:
            am_df['am_score'] = am_df['Score']
        else:
            am_df['am_score'] = am_df.iloc[:, 2]

        last_residue_position = am_df['Position'].max()

        # Cria o box plot dos scores preditos
        figAM2 = px.box(am_df, x='Position', y='am_score', title='Predicted AlphaMissense Scores per Position')
        figAM2.add_shape(type='rect', x0=0, x1=last_residue_position, y0=0, y1=0.33,
                         fillcolor='blue', opacity=0.1, line_width=0)
        figAM2.add_shape(type='rect', x0=0, x1=last_residue_position, y0=0.34, y1=0.56,
                         fillcolor='gray', opacity=0.1, line_width=0)
        figAM2.add_shape(type='rect', x0=0, x1=last_residue_position, y0=0.56, y1=1.0,
                         fillcolor='red', opacity=0.1, line_width=0)
        figAM2.show()
        save_graph(figAM2, "graph_AlphaMissense")
    else:
        print("O arquivo AlphaMissense não contém as colunas esperadas.")
else:
    print("Nenhum arquivo AlphaMissense carregado para análise.")


# In[ ]:





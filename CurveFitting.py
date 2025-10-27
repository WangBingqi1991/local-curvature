import numpy
import pandas as pd
import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import os
from itertools import cycle
import warnings

warnings.filterwarnings('ignore')

fusion_data = pd.read_excel(r'C:/Users/DAG/Desktop/bro/1025data/fusiondata.xlsx', header=None)
beauty_data = pd.read_excel(r'C:/Users/DAG/Desktop/bro/1025data/finaml1025image_all_stats.xlsx', header=None)

eta_true = pd.to_numeric(fusion_data.iloc[:, 0], errors='coerce').values / 100

if np.any(np.isnan(eta_true)):
    print("警告: fusion_data中存在非数值数据，已自动过滤")
    valid_indices = ~np.isnan(eta_true)
    eta_true = eta_true[valid_indices]

T = 2400

beauty_data.columns = ['k', 'i','d9595','numValidVesicles', 'flag', 'A_ik', 'perimeter'] + [f'hotpoint_{j}' for j in range(29)]

valid_data = beauty_data[beauty_data['flag'] == 1].copy()

numeric_columns = ['A_ik', 'perimeter'] + [f'hotpoint_{j}' for j in range(29)]
for col in numeric_columns:
    valid_data[col] = pd.to_numeric(valid_data[col], errors='coerce')

valid_data = valid_data.dropna(subset=numeric_columns)

unique_ks = valid_data['k'].unique()


def model_eta_single_k(lambda_val, p0, k_data, hotpoint_col, T_val=2400):

    numerator = 0
    denominator = 0

    for idx, row in k_data.iterrows():
        I_rho = row[hotpoint_col]
        A_ik = row['A_ik']

        if I_rho > 0 and A_ik > 0:
            exp_term = 1 - (np.exp(-I_rho * lambda_val * T_val)* (1 - p0))
            numerator += exp_term * A_ik
            denominator += A_ik

    if denominator == 0:
        return 0
    return numerator / denominator


def objective_function(params, all_k_data, eta_true, hotpoint_col):

    lambda_val, p0 = params
    predictions = []

    for k_data in all_k_data:
        eta_pred = model_eta_single_k(lambda_val, p0, k_data, hotpoint_col)
        predictions.append(eta_pred)

    predictions = np.array(predictions)

    min_len = min(len(predictions), len(eta_true))
    mse = np.mean((predictions[:min_len] - eta_true[:min_len]) ** 2)
    return mse


all_results = []
selected_k_data = []
k_ratios = []

for k in unique_ks:
    k_data = valid_data[valid_data['k'] == k]
    if len(k_data) > 0:
        selected_k_data.append(k_data)
        hotpoint_sum = k_data['hotpoint_10'].sum()
        perimeter_sum = k_data['perimeter'].sum()
        ratio = hotpoint_sum / perimeter_sum if perimeter_sum > 0 else 0
        k_ratios.append(ratio)
if len(selected_k_data) != 36:
    print('Wrong Column!')

for col_idx in range(6,29):
    hotpoint_col = f'hotpoint_{col_idx}'

    p0=numpy.random.random_sample()/100+1e-6

    initial_guess = [1e-6, p0]

    try:

        bounds = [(1e-10, 1e-2), (1e-6, 1e-2)]

        result = minimize(
            lambda params: objective_function(params,selected_k_data, eta_true, hotpoint_col),
            initial_guess,
            method='L-BFGS-B',
            bounds=bounds,
            options={'maxiter': 1000}
        )

        if result.success:
            lambda_opt, p0_opt = result.x

            predictions = []
            for k_data in selected_k_data:
                eta_pred = model_eta_single_k(lambda_opt, p0_opt, k_data, hotpoint_col)
                predictions.append(eta_pred)

            predictions = np.array(predictions)

            min_len = min(len(predictions), len(eta_true))
            predictions = predictions[:min_len]
            eta_used = eta_true[:min_len]

            ss_res = np.sum((eta_used - predictions) ** 2)
            ss_tot = np.sum((eta_used - np.mean(eta_used)) ** 2)
            r_squared = 1 - (ss_res / ss_tot) if ss_tot != 0 else 0


            all_results.append({
                'hotpoint_col': hotpoint_col,
                'lambda': lambda_opt,
                'p0': p0_opt,
                'r_squared': r_squared,
                'x_values': k_ratios,
                'predictions': predictions
            })

            print(f"热点列 {col_idx+1}: λ={lambda_opt:.6e}, p0={p0_opt:.6f}, R²={r_squared:.6f}\n，预测值是：")
            for i in predictions:
                print(f"{i:.6f}",end=' ')
            print()
        else:
            print(f"热点列 {col_idx+1} 优化失败: {result.message}")

    except Exception as e:
        print(f"热点列 {col_idx+1} 拟合失败: {e}")
        continue

if all_results:
    all_results.sort(key=lambda x: x['r_squared'], reverse=True)
    best_results = all_results[:5]

    save_dir = r'C:\Users\DAG\Desktop\图片'
    os.makedirs(save_dir, exist_ok=True)

    colors = cycle(['blue', 'red', 'green', 'orange', 'purple'])

    print("\n=== 最好的五个拟合结果 ===")

    for i, result in enumerate(best_results):
        color = next(colors)

        print(f"\n结果 {i + 1}:")
        print(f"热点列: {result['hotpoint_col']}")
        print(f"λ = {result['lambda']:.6e}")
        print(f"p0 = {result['p0']:.6f}")
        print(f"R² = {result['r_squared']:.6f}")

        plt.figure(figsize=(10, 6))

        plt.scatter(result['x_values'], result['eta_true'], color='black', label='真实值', alpha=0.7, s=50)

        plt.scatter(result['x_values'], result['predictions'], color=color, label='预测值', alpha=0.7, s=50)

        plt.plot(result['x_values'], result['eta_true'], 'k--', alpha=0.5, linewidth=1, label='真实值趋势')
        plt.plot(result['x_values'], result['predictions'], color=color, linestyle='-', linewidth=2, label='预测值趋势')

        plt.xlabel('hotPointNums / perimeter')
        plt.ylabel('η')
        plt.title(f'拟合结果 {i + 1} (R² = {result["r_squared"]:.4f})')
        plt.legend()
        plt.grid(True, alpha=0.3)

        plt.savefig(os.path.join(save_dir, f'best_fit_{i + 1}.png'), dpi=300, bbox_inches='tight')
        plt.close()

    print(f"\n图片已保存到: {save_dir}")

    print("\n=== 最终结果汇总 ===")
    for i, result in enumerate(best_results):
        print(f"\n结果 {i + 1}:")
        print(f"λ = {result['lambda']:.6e}")
        print(f"p0 = {result['p0']:.6f}")
        print(f"R² = {result['r_squared']:.6f}")
else:
    print("没有成功的拟合结果，请检查数据格式和模型参数")
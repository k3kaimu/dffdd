from ast import dump
import math
import optuna
import subprocess
from multiprocessing import Process


def objective(trial):
    # """最小化する目的関数"""
    # # パラメータが取りうる範囲
    # x = trial.suggest_uniform('x', -5, +15)
    # # デフォルトで最小化かつ現在は最小化のみのサポートなので符号を反転する
    # return - math.exp(-(x - 2) ** 2) + math.exp(-(x - 6) ** 2 / 10) + 1 / (x ** 2 + 1)
    # inner_LDPC_iter = trial.suggest_int('inner_LDPC_iter', 1, 5)
    # inner_GaBP_iter = trial.suggest_int('inner_GaBP_iter', 1, 5)
    damp1_1 = trial.suggest_uniform('damp1_1', 0, 1)
    damp1_2 = trial.suggest_uniform('damp1_2', 0, 1)
    damp2_1 = trial.suggest_uniform('damp2_1', 0, 1)
    damp2_2 = trial.suggest_uniform('damp2_2', 0, 1)
    # numiter = trial.suggest_int('numiter', 1, 10)

    ret = subprocess.run(["./gabp_gaussian", "50", str(damp1_1), str(damp1_2), str(damp2_1), str(damp2_2)], capture_output=True)
    # ret = subprocess.run(["./gabp_gaussian", str(damp1_1)], capture_output=True)
    return float(ret.stdout)



def optimize(study_name, storage, n_trials):
    # 最適化のセッションを作る
    sampler = optuna.samplers.CmaEsSampler()
    # sampler = optuna.samplers.RandomSampler()
    study = optuna.load_study(
        study_name=study_name,
        storage=storage,
        sampler=sampler
    )
    # study.enqueue_trial({'damp1_1': 1, 'damp1_2': 0.05, 'damp2_1': 0.2, 'damp2_2': 0.2})
    study.enqueue_trial({'numiter': 20, 'damp1_1': 0.15853341118873465, 'damp1_2': 0.023387205571840916, 'damp2_1': 0.5915693314590689, 'damp2_2': 0.6479578055604401})
    study.optimize(objective, n_trials=n_trials)


if __name__ == '__main__':
    DATABASE_URI =  'sqlite:///optuna_db.db'
    study_name = 'sefdm_gabp_ldpc'

    study = optuna.create_study(storage=DATABASE_URI, study_name=study_name)

    workers = [
        Process(
            target=optimize, 
            args=(study_name, DATABASE_URI, 1000)
        ) for _ in range(1)
    ]

    for worker in workers:
        worker.start()

    for worker in workers:
        worker.join()

    # 最適化の結果を確認
    study = optuna.load_study(study_name=study_name, storage=DATABASE_URI)
    best_trial = study.best_trial
    print(f"""Number of finished trials: {len(study.trials)}
    Best trial:
        Value: {best_trial.value:.6f}
    Params: """)
    for k, v in best_trial.params.items():
        if type(v) is float:
            print(f'        {k}: {v:.6f}')
        else:
            print(f'        {k}: {v}')

    # optuna.visualization.plot_contour(study)
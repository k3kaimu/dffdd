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
    inner_LDPC_iter = trial.suggest_int('inner_LDPC_iter', 1, 5)
    inner_GaBP_iter = trial.suggest_int('inner_GaBP_iter', 1, 5)
    dumping = trial.suggest_uniform('dumping', 0, 1)
    outer_iter = trial.suggest_int('outer_iter', 10, 30)

    ret = subprocess.run(["./gabp_gaussian", str(inner_LDPC_iter), str(inner_GaBP_iter), str(dumping), str(outer_iter)], capture_output=True)
    return float(ret.stdout)



def optimize(study_name, storage, n_trials):
    # 最適化のセッションを作る
    # sampler = optuna.samplers.CmaEsSampler()
    sampler = optuna.samplers.RandomSampler()
    study = optuna.load_study(
        study_name=study_name,
        storage=storage,
        sampler=sampler
    )
    study.optimize(objective, n_trials=n_trials)


if __name__ == '__main__':
    DATABASE_URI =  'sqlite:///optuna_db.db'
    study_name = 'sefdm_gabp_ldpc'

    workers = [
        Process(
            target=optimize, 
            args=(study_name, DATABASE_URI, 1000)
        ) for _ in range(10)
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
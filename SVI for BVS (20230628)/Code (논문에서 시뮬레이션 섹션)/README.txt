Code (for simulation section of our paper)
- 논문 SVI for BVS (20230628)의 Simulation section과 Real data application section에 쓰인 코드들

1. 4개의 R 파일들 중 가장 먼저 숙지할 코드:
- synthetic_data_analysis_high_dim.R
- 왜냐면 여기에 중요한 코드들에 대한 주석들이 있기 때문
- 이 코드들을 익히면 나머지 synthetic_data_analysis_low_dim.R과 real_data_analysis.R들도 수월히 이해할 수 있음

2. synthetic_data_analysis_high_dim.R의 전반적 흐름 :
- 라이브러리 / 전역함수 / Synthetic data 생성함수 / 알고리즘 / 시뮬레이션세팅 / 시뮬레이션실행

3. 주의사항 :
- 시뮬레이션이 있는 3개의 R 파일들 모두 MFVI 알고리즘 작성에서 stopping criterion (pred_hist) 코드 수정을 안함. ELBO로 했어야하나 현재는 RMSE로 한 상황
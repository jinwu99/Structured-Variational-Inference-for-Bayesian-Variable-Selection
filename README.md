### 0. 지금까지 진행상황

* SVI for BVS (20230209)
* SVI for BVS (20230628)



### 1. Variable selection in Regression Problem

$$
y=\sum_{j=1}^p x_j\beta_j + \epsilon, \quad \epsilon \sim N(0,\sigma^2)
$$

* 대표적인 변수선택 방법 3가지 : Criteria-based methods(AIC/BIC...), Penalized Likelihood(Lasso/SCAD...), Bayesian approach(shrinkage prior/spike and slab prior)

* Spike and slab prior :
  $$
  \beta_j|\gamma_j,\sigma^2 \sim \gamma_j \pi(\theta_j|\sigma^2) + (1-\gamma_j)\delta_0 (\beta_j) \\
  \gamma_j = \begin{cases} 1 & \text{if variable $j$ is included in model,} \\
  0 & \text{o.w.} \end{cases}
  $$

  * Slab distribution $\pi(\theta_j|\sigma^2)$ : 변수가 포함될 때$(\beta_j \neq 0)$의 분포.
    * Noise variance  $\sigma^2$도 prior distribution에 포함하는 이유 : 데이터의 노이즈 강도가 회귀계수 추정에도 영향미칠것을 고려
  * Spike distribution $\delta_0$ : $\beta_j=0$에서 Dirac function. Sparsity induce해줌
  * Gold standard for high-dimensional(변수가 데이터보다 많은 경우) variable selection
  * 직관적인 해석 : $E[\gamma_j|y,x_1,...,x_p]$ : 데이터가 주어졌을 때 $x_j$가 유효할 사후확률. 주로 1/2보다 크면 변수선택함

* Which distribution to use for slab prior?

  * Gaussian slab $N(\theta_j|0,\sigma^2\sigma_{\theta_j}^2)$ : 사후분포에 대한 계산이 편리해서 많이 쓰이지만 excessive shrinkage를 일으킴, 이에 따라 변수선택 및 MSE에 악영향 가능.
  * Laplace slab $Laplace(\theta_j|0,\lambda/\sigma)$ : Gaussian slab보다 두꺼운 꼬리를 가져서(exponential tail) 회귀계수가 extreme value인 경우에도 보다 많은 mass를 부여할 수 있음. 따라서 사후분포가 더 빨리 수렴하고 MSE 및 F1-score에도 긍정적인 영향 미칠 것으로 기대. 하지만 지나치게 큰 $\lambda$값은 shrinkage를 지나치게 일으킬 위험 있어 적당한 값을 줘야함. 현재 연구에서는 이 분포를 사용하며 $\lambda=1$로 fix하여 사용하는중.



#### 1.1. Model Likelihood and Prior distributions

$y=(y_1,...,y_n)^\top$와 $n\times p$사이즈의 design matrix $X$에 대하여 
$$
\begin{align}
p(y|\gamma,\theta,\sigma^2) &= N(y|X\Gamma\theta,\sigma^2 I_n), \quad \pi(\sigma^2) \propto 1/\sigma^2, \\
\pi(\theta| \tau^2,\sigma^2) &= N(\theta|0,\sigma^2\text{diag}(\tau^2)), \quad \pi(\tau^2) = \prod_{j=1}^p Exp(\tau^2_j|\lambda^2/2) \\
\pi(\gamma|\rho) &= \prod_{j=1}^p Bernoulli(\gamma_j|\rho), \quad \pi(\rho) = Beta(A_0,B_0)
\end{align}
$$

* $y=(y_1,...,y_n)^T$ : continuous data, 

  $X : n\times p$ design matrix,

  $\theta = (\theta_1,...,\theta_p)^T$ : coeffieicnts

  $\gamma=(\gamma_1,...,\gamma_p)^T$ : binary indicators for variable selection, $\Gamma = \text{diag}(\gamma)$

  $\tau^2=(\tau^2_1,...,\tau^2_p)^T$ : variance of $\theta$,

  $Exp(\lambda^2/2)$ : exponential distribution

* Laplace slab and spike prior : $\theta$의 prior 분포에 대하여 각 $\tau_j^2$에 marginalize out하면 각 변수 $\beta_j=\gamma_j\theta_j$는 laplace slab and spike prior을 갖는다 :
  $$
  \pi(\theta_j|\sigma^2,\lambda) = \int N(0,\sigma^2\tau_j^2)Exp(\lambda^2/2) d\tau^2_j = Laplace(\theta_j|0,\lambda/\sigma),
  \\ \beta_j|\gamma_j \sim \gamma_j Laplace(\theta_j|0,\lambda/\sigma) + (1-\gamma_j)\delta_0
  $$

* $\sigma^2$에 improper prior : noninformativeness를 위해 부여, 모든 가능한 모델$\mathcal{M}_\gamma$에 대해 공통적으로 들어가있는 파라미터라서 improper prior부여 가능.

* 변수선택 결과의 sparsity에 아주 큰 영향미치는 파라미터가 laplace slab의 $\lambda$와 배르누이분포의 $\rho$임. $p>>n$인 high-dimensional variable selection에서 많은 경우 $\pi(\rho)$의 셋팅을 $Beta(1,p)$로 함(이론적으로도 실용적으로도 좋다고 알려져있음) 따라서 그렇게 설정함.



#### References

Bayesian Variable Selection 폴더 안에 README 참조.



### 2. Variational Inference (VI)

대략적 이론 설명 :

* 잠재변수 $z$의 사후분포를 직접 구하기 어려운 경우가 많음(due to normalizing constant, etc). MCMC로 사후분포의 샘플들을 얻을 수 있지만 고차원 데이터에 대해서는 수렴이 오래걸리고 수렴성 진단도 쉽지 않음. VI는 사후분포를 가장 최적으로 근사하는 변분분포를 구하는 방법인데
  $$
  q^*(z) = \arg\min_{q\in Q} D_{KL}(q(z)\|p(z|y)), \\
  \text{where }D_{KL}(q(z)\|p(z|y)) = E_{q(z)}\left[\log \frac{q(z)}{p(z|y)}\right]
  $$
  사후분포 $p(z|y)$와의 Kullback leibler(KL) divergence를 최소화하는 최적의 변분분포 $q(z) \in Q$를 구하는 방법. 그런데 직접 계산하기 어려움으로 주로 Evidence lower bound(ELBO)를 최대화하는 방법으로 구함.
  $$
  \begin{align}
  q^*(z) &= \arg\min_{q\in Q} D_{KL}(q(z)\|p(z|y)), \\
  &= \arg \max_{q\in Q}\mathcal{L}(q), \\
  \text{where }&\mathcal{L}(q)=\underbrace{E_{q(z)}[\log p(z,y)] + H(q)}_{ELBO} \\
  \text{and }&H(q) = \underbrace{E_{q(z)}[-\log q(z)]}_{entropy}
  \end{align}
  $$
  
* 분포족 $Q$을 어떤 걸로 정하느냐 따라 VI알고리즘의 난이도가 달라짐.

  * 보통 독립가정을 많이 할수록 최적분포를 tractable하게 구하기 쉬워지되 근사정확도가 떨어짐. 반면 독립가정을 덜 할수록 근사정확도가 비교적 좋아지지만 최적분포를 tractable하게 구하기 어려울 수 있음.

    (애초에 사후분포를 exact하게 구하기 어려운 이유가 잠재변수들간의 dependence임을 생각해보면 일부러 독립가정을 넣어 근사를 쉽게 하려는 목적임을 알 수 있음)

* 만약 분포족을 $Q=\{q:q(z)=\prod_{l=1}^q q(z_l)\}$로 한정시킨다면 (이때 각 $z_l$은 벡터일 수 있음) 각 최적의 분포$q^*(z_l)$를 구하는 방법은 다음과 같은 공식으로 구할 수 있음 :
  $$
  q^*(z_l) \propto E_{q(z_{-l})}[\log q(z_l|z_{-l})]
  $$

  * 이때 $z_{-l}$은 $z_l$을 제외한 모든 잠재변수.
  * Gibbs sampling의 full conditional form과 비슷함(단지 log씌우고 $q(z_{-l})$에 대한 기대값 씌워준거에 차이가 있음) 따라서 Gibbs sampling처럼 각 잠재변수들의 prior과 liklihood간에 전부 exponential family이고 conjugacy 만족하면 쉽게 form을 구할 수 있음.
  * 업데이트과정은 coordinate ascent algorithm : 목적함수인 범함수 ELBO에 대하여 각 coordinate $q(z_l)$을 순차적으로 업데이트. 
    * ex. $q(z_1,z_2,z_3)=q(z_1)q(z_2)q(z_3)$이라면 $q(z_1)\rightarrow q(z_2)\rightarrow q(z_3)$. 이때 업데이트 순서는 달라질 수 있음.

* 모든 잠재변수에 대해 전부 다 독립으로 가정하면 Mean Field VI(MFVI)라고 한다. 우리 모형의 경우 MFVI모형의 가정을 다음과 같이 정의함 :
  $$
  q(\theta,\gamma,\sigma^2,\rho,\tau^2) = q(\theta)\left\{\prod_{j=1}^p q(\gamma_j)\right\}q(\sigma^2)q(\rho)\left\{\prod_{j=1}^p q(\tau_j^2)\right\}
  $$

  * 각 최적분포의 form과 유도과정은 다음 자료를 참고바람 : SVI for BVS (20230209)



#### 2.1. Determinstic Annealed VI (DAVI)

대략적 이론 설명 :

* VI의 문제점은 주로 목적함수 ELBO가 nonconvex하다는 것임(starting point와 변분분포들의 업데이트 순서에 따라 local optima에 갇힐 염려있음). 특히 독립가정을 많이 할수록 nonconvexity가 강해짐.

* 이를 해결하는 방법 중 하나가 Deterministic Annealing을 활용하는 것. 맨 처음에 ELBO를 annealing을 통해 강제로 concave하게 만들고 업데이트를 진행함에 따라 annealing을 줄여주어 원래의 ELBO로 돌아가게 하는것

* 수식을 통해 본다면 ELBO는 변분분포들 뿐만 아니라 annealing parameter $T$에 대해서도 의지한다 :
  $$
  \mathcal{L}(q;T) = E_{q(z)}[\log p(z,y)] + T\times H(q)
  $$

  * 엔트로피 $H(q)$는 분포 $q(z)$에 대하여 concave하다. 온도와 같은 역할 하는 파라미터$T$를 크게 할수록 $\mathcal{L}(q;T)$ 또한 concave형태에 가깝다. 그러면 coordinate ascent 과정에서 $\mathcal{L}(q;T)$의 global optima에 갈 가능성 또한 커진다.

  * 온도 $T$를 점차 내려주며 업데이트 진행하다가 $T=1$이 되었을때 starting point는 보다 좋은 local optima로 가기 위한 상태일 것이다.

    <img src="C:\Users\이진우\AppData\Roaming\Typora\typora-user-images\image-20230807155055059.png" alt="image-20230807155055059" style="zoom:60%;" />

  Variational Tempering 논문에서 가져온 그림.
  
* 우리 알고리즘의 경우 $q(\gamma)=\prod_{j=1}^p q(\gamma_j)$에 대해서만 annealing을 주었다 :
  $$
  \mathcal{L}(q;T) = E_{q(z)}[\log p(z,y)] + H(q(z_{-\gamma})) + T\times H(q(\gamma))
  $$
  
  * 이유는 잘 모르겠으나 다른 변분분포들에 대해서도 함께 annealing을 주었을 때는 잘 작동하지 않음. (항상 Null모형 - 아무 변수도 선택안한 모형만 최적분포로 도출됨)



#### References

VI Variants 폴더 안에 README 참조.



#### 2.2. Structured VI (SVI)

대략적 이론 설명 :

* MFVI의 완전독립가정을 완화한 것이 SVI. 통계적 의존성을 회복하는 만큼 실제 사후분포에 근사를 더 잘할 것으로 기대. ELBO도 nonconvexity가 완화되어 initialization에 덜 민감해지는 것으로 실험적으로 알려져있음.
  * 다만 의존성을 극복하는 데에 tractable한 최적분포가 도출되지 못하는 cost가 있을 수 있음.

* **References** : VI Variants 폴더 안에 README 참조.
  
* 지금까지 두 가지 Structured Variational form을 제시함 :

$$
\begin{align}
&(a) : q(\theta,\gamma,\sigma^2,\rho,\tau^2) = q(\theta)q(\gamma)q(\sigma^2)q(\rho)\left\{\prod_{j=1}^p q(\tau_j^2)\right\} \\
&(b) : q(\theta,\gamma,\sigma^2,\rho,\tau^2) = q(\theta,\gamma,\sigma^2)q(\rho)\left\{\prod_{j=1}^p q(\tau^2_j)\right\}
\end{align}
$$



##### 2.2.(a).

Derivation 참고 : SVI for BVS (20230209)

* $q(\gamma)$은 Ising 모형과 비슷한 form을 가진다 :
  $$
  q(\gamma) = \frac{Q(\gamma)}{Z},\quad Q(z)=\exp(\psi^T\gamma + \gamma^T\Psi \gamma), \quad Z=\sum_{\gamma}Q(z)
  $$

  * Graph에서 $\gamma$를 node들이라 볼 때 $\psi$는 node들의 potential, $\Psi$는 node들간의 edge strength(Adjacency matrix)이다. 
  * MFVI가정과 달리 $\gamma$노드들간에 dependency가 고려된다.

* 문제점 : 변분분포들을 업데이트하기 위해 $E_{q(\gamma)}[f(\gamma)]$계산이 필요한데 이는 $2^p$만큼의 summation을 요구로한다. 즉, $p$가 작을 때는 exact computation이 가능하지만 커질수록 불가능해짐. 기대값을 approximate하는 방법을 생각해야함.



##### Ideas for approximating $E_{q(\gamma)}[f(\gamma)]$:

##### 1. Approximate with Gibbs Sampling (GS)

$$
E_{q(\gamma)}[f(\gamma)] \approx \frac{1}{S}\sum_{s=1}^S f(\gamma^{s}), \quad \{\gamma^s\}_{s=1}^S \sim q(\gamma)
$$

장점 :

1. $2^p$만큼 연산할 필요가 없음. GS를 통해 몇몇 high-probability GS samples $\{\gamma^s\}$만 활용하여 기댓값 근사 가능

단점 :

1. 변수들간의 correlation이 커질수록 $q(\gamma)$에 대한 GS의 mixing이 잘 안됨. GS가 잘 수렴하지 않은 상태에서의 $\{\gamma^s\}$들을 사용하면 오히려 approximate 정확도가 떨어지는데, 심지어 MFVI보다도 정확도가 안좋아짐.

2. 변수 갯수가 많아질수록 GS mixing에 필요한 GS횟수가 더 커질 수 있음, 속도저하 (VI는 빠른 속도가 매우 중요!)



##### 2. Approximate with Sequential Monte Carlo Sampling (SMCS) :

기댓값을 단순 GS대신 Sequential Monte Carlo Sampler(with multiple GS chains)으로 근사

SMCS의 대략적 설명 :

* MCMC는 주로 chain하나만 활용하여 샘플링하는데, 특히 GS는 local mode에 갇힐 위험성이 있음.
* SMCS는 여러 개의 chain들을 활용하고 중간중간 Importance Sampling(IS)을 통해 여러 chain의 샘플들 중 중 좋은 샘플만 추려내는 과정이 있다.
* 자세한 설명은 SVI for BVS (20230209)와 밑의 참고자료들 보기 바람.
  * **References** : VI Variants 폴더 안에 README 참조.


장점 :

1. GS의 mixing문제를 어느정도 해결할 수 있을 것으로 기대함.

단점 :

1. 근데 생각보다 GS보다 성능이 좋게 나오는 경우가 많지 않고 / 심지어 좋게 나와도 성능향상은 미비하고 시간은 GS보다도 느림 (병렬 컴퓨팅을 활용했음에도 그러함.)

*개인적으로 시간을 많이 쏟은 방법이지만 효과는 크지 않기에 SMCS자체 방법에 흥미를 느끼시지 않는 이상 그냥 넘어가시는걸 추천합니다...*



##### 3. Approximate with Monte Carlo sampling - but with no MCMC sampling :

MCMC의 그 어떤 방법을 활용하든 Large discrete space안에서 효율적으로 mixing되는건 지금도 아주아주 어려운 문제임.

따라서 차라리 MCMC처럼 dependent sampling이 아닌 순수 indepdent sampling이 가능한 단순 샘플러를 활용하고 그 샘플러들을 점차 최적화하는 방법도 생각해 볼 수 있음.

다양한 방법들 고안해낼 수 있는데 - 지금까지 나온 아이디어들은 :

1. VI 알고리즘에서 $i$번째 업데이트된 $Q_{[i]}(\gamma)$의 기댓값을 계산하기 위하여 :

   * 처음에 Uniform distribution에서 엄청나게 많이 샘플링 (e.g.  $\{\gamma^s_{[i]}\}_{s=1}^S\sim Bernoulli(1/2)$)을 한 이후 $Q_{[i]}(\gamma)$를 활용하여 Importance sampling(IS)을 한다.
   * IS르르 통해 필터되고 남은 unique 샘플 종류가 극히 적을 것이기에 이 샘플들에 조금씩 노이즈를 섞고 다시 IS한다. 여러번 반복 (이때 매번 IS할때마다 필터링된 샘플들은 계속 저장해나간다)
   * 지금까지 저장된 샘플들을 활용해 기댓값 $E_{q_{[i]}(\gamma)}[f(\gamma)]$ approximate

2. VI 알고리즘에서 $i$번째 업데이트된 $Q_{[i]}(\gamma)$의 기댓값을 계산하기 위하여 :

   * 처음에 Uniform distribution에서 엄청나게 많이 샘플링 (e.g.  $\{\gamma^s_{[i]}\}_{s=1}^S\sim Bernoulli(1/2)$)을 한 이후 $Q_{[i]}(\gamma)$를 활용하여 IS를 한다.

   * 필터되고 남은 샘플들을 활용하여 샘플러 $Bernoulli(1/2)$갱신 (e.g. $Bernoulli(\hat\phi^{mle})$).

   * 갱신된 샘플러를 활용해 또 엄청나게 많이 샘플링하고 $Q_{[i]}(\gamma)$활용하여 IS다시한다. (몇번 반복할지는 선택사항)

   * 지금까지 저장된 샘플들을 활용해 기댓값 $E_{q_{[i]}(\gamma)}[f(\gamma)]$ approximate

그러나 이 두가지 모두 뚜렷한 성과가 없던 것으로 기억...



##### 4. 굳이 $\gamma$를 전부 structured하지 말고 일부 혹은 partition별로 하자. 

예를 들어
$$
q(\gamma) = q(\gamma_{1:10})\prod_{j=11}^p q(\gamma_j)
$$
굳이 순차적인 인덱스 $\gamma_{1:10}$로 묶는게 아니라 다양한 criterion으로 묶을 수 있다. (예를 들어 Correlation이 가장 높은 10개의 변수들끼리 묶기)

또는
$$
q(\gamma) = q(\gamma_{1:10})q(\gamma_{11:20})q(\gamma_{21:30})...
$$
이것 또한 Correlation순위대로 partition시킬 수 있다.

장점 :

1. GS효율성이 커진다 : 이전처럼 p-차원의 discrete space에서 MCMC가 mixing하는건 너무 힘든데 그에 비해 10-차원에서 mixing하는건 훨씬 수월하기 때문.

단점 :

1. 어떻게 묶느냐가 문제일 수 있다. 예로 Correlation이 높은 순서대로 묶으면 Correlation이 높은 변수들이 함께 선택될 가능성이 높다 - 유전자 데이터처럼 연관된 유전자가 함께 선택되어야할 경우 유리하나 Multicollinearity가 심한 경우에도 중요한 변수만 선택하고 싶은 경우엔 불리함.
2. MFVI의 문제점을 여전히 가지고 있음 : 어떤 partition $q(\gamma_{j:{j+c}})$을 먼저 업데이트하냐에 따라 local optima가 달라짐 - VI는 독립가정이 강할수록 ELBO의 nonconvexity가 강해지기 때문

*뚜렷한 성과가 없던 것으로 기억나는데 제가 전에 만들었던 코드가 엉망이 되어서ㅠ 시도해보신다면 직접 코드 만드시길 부탁드립니다.*



##### 2.2.(b).

Derivation 참고 : SVI for BVS (20230628)

* $q(\gamma)$ 분자의 $Q(\gamma)$연산에 matrix inverse와 matrix determinant계산이 엄청 빡쎄다. 즉 **2.2.(a).**에서 논의했던 approximation 방법들 **2,3**번 모두 힘듬, 왜냐면 $Q(\gamma)$연산을 고반복해야하기 때문이다. 
* 4번의 approximation방법 경우 각 $q(\gamma_{j:j+c})$를 유도하기 위해 expectation of matrix inverse와 expectation of matrix determinant 계산이 가능해야하는데 아마 불가능할 것임. 따라서 4번도 불가능.
* 그래서 지금까지 쓴 논문에는 1번 approximation(GS approximate)만 활용하고 있음.



### 3. 당부사항

* 저는 variable selection에 VI를 적용하는 문제를 위주로 다루었지만 다음 연구원분께서는 새로운 VI 응용방법을 만들거나 최근 VI방법을 활용하여 다른 문제를 해결하는 방향으로 가시는 것도 좋다고 생각합니다. VI variant폴더를 통해서도 여러 VI응용방법들을 보실 수 있습니다.

* 그래도 만약 미완성 논문 SVI for BVS (20230628)의 연구를 이어나가시게 된다면 먼저 SVI for BVS (20230628)에 해당하는 R코드들부터 확인해주시기 부탁드립니다. 해당 코드들이 가장 정제되고 주석이 많습니다.

  숙지되신 이후 만약 $2.2.(a)$의 여러 approximation 방법에도 관심있으시다면 그때 SVI for BVS (20230209) 코드를 확인해주세요.

* 질문이 있으시면 제게 카톡이나 메세지로 연락해주시기 바랍니다 : 010-7579-5686, 이진우

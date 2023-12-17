# Introduction
The Expectation Maximization (EM) algorithm is a method used for estimating the parameters of a statistical model that includes hidden variables. In each iteration, the EM algorithm consists of two steps: the Expectation step (E-step) and the Maximization step (M-step). The E-step calculates the expected value of the log-likelihood given the current estimate of the model parameters, while the M-step updates the parameters to maximize this expected log-likelihood. This is similar to being in a valley where the whole mountain can’t be seen, but its shape can be estimated based on current observations, which determined the next step to go. This is somewhat similar to the gradient descent algorithm in machine learning.

<div align="center">
  <img src="https://github.com/fumiama/protein-motifs-em-algorithm/assets/41315874/c5073503-2ff2-42a4-9d6f-ea74ac8f1ac7">
</div>

# Usage
Just download the whole repo and run `em.m` to see the result.

# Thanks
- [Furuzuki Laboratory, IPS, Waseda University](https://nclab.w.waseda.jp/nclab/index.html)

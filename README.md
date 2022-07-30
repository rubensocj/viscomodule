# viscomodule

The *prev-NAS* is a R project for modeling and forecasting groundwater level (GWL) time series with combined models.

*Read this in other languages*: [Português brasileiro](README.br.md).

## Theorical background

A partir da integral de convolução para a tensão no instante $t$, da série de Prony para o Módulo de Relaxação e da expressão de entrada do ensaio com taxa de deformação constante $k_0$ com o tempo $\epsilon(t) = k_0 t$, a função de interpolação utilizada é dada a seguir.

$$\sigma (t) = k_0 \big[ E_{\infty t} + \sum_{i=1}^{n} E_i \rho_i (1 - e^\frac{-t}{\rho_i}) \big]$$

em que $n$ é o número de termos da série de Prony e os valores de $\sigma$ são obtidos do ensaio com taxa de deformação constante.

## Installation

You can download this source code and simply run the [viscomodule](viscomodule.py) file.

## Example

Relaxation module from the [literature](relaxation-module) and real data from [CCMRT test](CCMRT-test/CCMRT-test.csv) at UFCA are avaliable.

## Citation
If you use the code of this repository in your paper or research please cite:

```
@MASTERSTHESIS{,
  title = {Determinação das constantes da série de Prony do módulo de relaxação de misturas asfálticas por ensaio de compressão uniaxial com taxa de deformação constante}
  author = {da Cunha Junior, Rubens Oliveira},
  year = {2019},
  school = {Universidade Federal do Cariri},
  address = {Juazeiro do Norte},
  type = {Graduação em Engenharia Civil}
}
```

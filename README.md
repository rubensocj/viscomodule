# viscomodule

O *viscomodule* é um projeto em `Python` para caracterização viscoelástica do módulo de relaxação de misturas asfálticas usando séries de Prony.

A partir dos resultados de um ensaio de _creep_ estático, o código considera apenas o regime viscoelástico linear e permite obter as constantes da série de Prony do módulo de relaxação. Uma função linear é usada como função interpoladora para realizar o ajuste de curvas da série de Prony e o método dos mínimo quadrados é aplicado para se obter o sistema de equações lineares.

Usando uma série de Prony do módulo de relaxação previamente conhecida, o programa também permite simular um ensaio de _creep_ estático.

Módulos de relaxação da [literatura](relaxation-modulus) e dados reais de um ensaio de [_creep_ estático](creep-test/creep-test.csv) realizado na UFCA estão disponíveis.

O projeto possui uma Interface Gráfica do Usuário (GUI) em `PyQt5`.

Este projeto foi desenvolvido para o Trabalho de Conclusão de Curso (TCC) em Engenharia Civil na Universidade Federal do Cariri (UFCA).

## Instalação

Faça download do código e execute o arquivo [viscomodule](viscomodule.py).

## Citação

```
@MASTERSTHESIS{,
  title = {Determinação das constantes da série de Prony do módulo de relaxação de misturas asfálticas por ensaio de compressão uniaxial com taxa de deformação constante}
  author = {da Cunha Júnior, Rubens Oliveira},
  year = {2019},
  school = {Universidade Federal do Cariri},
  address = {Juazeiro do Norte},
  type = {Graduação em Engenharia Civil}
}
```

# Desborramento de imagens
Códigos desenvolvidos durante a disciplina **PGCC004 - Problemas Inversos Discretos (Teoria e Aplicação)**.
## Resumo
 O repositório se concentra em reproduzir técnica de desborramento de imagens utilizando regularização de Tikhonov, utilizando abordagem proposta por (Kamm, Nagy, 1998) [1]. A técnica é testada utilizando o conjunto de dados "I Can See Clearly Now" [2], utilizando GCV para escolha do parâmetro de regularização. Por fim, a abordagem de (Kamm, Nagy, 1998) é comparada com o método dos gradientes conjugados (CGLS).

## O problema 
Uma imagem borrada $\mathbf{G}$ pode ser representada pela convolução de uma *Point Spread Function* (PSF) com a imagem nítida $\mathbf{F}$, que deseja-se encontrar
$$
\mathbf{G}=\mathbf{P}*\mathbf{F}.
$$

Em forma de sistema linear, as colunas de $\mathbf{G}$ e $\mathbf{F}$ são empilhadas em um único vetor coluna por meio do operador $\mathrm{vec}$, e $\mathbf{P}$ se torna um tipo de matriz conhecido como BTTB(Bloco Toeplitz com Blocos Toeplitz), que cresce quadráticamente com as dimensões de $\mathbf{P}$, resultando na equação
$$
\mathbf{g} = \mathbf{K}\mathbf{f},
$$
em que:
- $\mathbf{g} = \mathrm{vec}(\mathbf{G})$
- $\mathbf{f} = \mathrm{vec}(\mathbf{F})$
- $\mathbf{K}_{n\times n} = \mathrm{toep2}(\mathbf{P}_{\sqrt{n}\times\sqrt{n}})$

Haja visto que $\mathbf{G}$ é enviesado por ruído e $\mathbf{K}$ apresenta grandes dimensões (considerando uma PSF de 256x256 pixels float64, a matriz ocupará 32 GiB), o problema é considerado **mal-posto**, criando a necessidade de técnicas de regularização para o problema.

(Kamm, Nagy, 1998) propõe representar $\mathbf{K}$ pelo produto de Kronecker de duas matrizes $\mathbf{A}$ e $\mathbf{B}$, evitando alocação direta de $\mathbf{K}$
$$
\mathbf{K} = (\mathbf{A}\otimes\mathbf{B}),
$$
dessa forma, a SVD de $\mathbf{K}$ pode ser representada em termos da SVD de $\mathbf{A}$ e $\mathbf{B}$
$$
\mathbf{K} = \mathbf{U}\Sigma \mathbf{V}^\mathrm{T} = (\mathbf{U}_A \otimes \mathbf{U}_B) (\Sigma_A \otimes \Sigma_B) (\mathbf{V}_A \otimes \mathbf{V}_B)^\mathrm{T},
$$
permitindo assim o uso de técnicas como regularização de Tikhonov e TSVD.


## Coisas que faltam
- s3_TikhonovDeducao
  - Plot de picard meio iffy
  - Como que eu chego no produto de hadamard

## Referências
[1] KAMM, J.; NAGY, J. G. Kronecker product and SVD approximations in image restoration. Linear Algebra and its Applications, v. 284, n. 1–3, p. 177–192, 1998. DOI: 10.1016/S0024-3795(98)10024-1.

[2] NAGY, J. G.; O’LEARY, D. P. Image deblurring: I can see clearly now. Computing in Science & Engineering, v. 5, n. 3, p. 82–84, 2003. DOI: 10.1109/MCISE.2003.1196312.

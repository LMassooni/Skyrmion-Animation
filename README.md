# Skyrmion-Animation
Animação de um sistema de spins sujeitos a interação de troca e de Dziyaloshinskii-Moryia
Basicamente, usei o gnuplot para fazer uma animação em tempo real da evolução de um sistema sujeito a interação de troca ferromagnetica e a interação anti-simétrica (DM) 
usando simulações de Monte Carlo, com Algoritmo de Metrópolis. No código é usado também um algorítmo de over-relaxation para acelerar a dinâmica das estruturas moduladas que
emergem da competição de ambas as interações: Um skyrmion.

A depender dos parâmetros internos (temperatura T, campo magnético externo H e anisotropia B) é possível termos diferentes fases moduladas. A primeira chama-se fase helicoidal
caracterizada por faixas em duas dimensões. A segunda é a fase de skyrmions, uma fase de "bolhas" caracterizada por um centro antiparalelo com as bordas, com um padrão de redemoinho entre
ambas.

O arquivo é auto-consistente. Basta baixá-lo, compilar com:

gcc animacao.c -lm -O3

e executar o executável, passando para o gnuplot com um pipe:

./a.out | gnuplot

Modificações nos plots podem ser feitos em duas funções dentro do código:
initGnuplot() -> Possui as características do Plot: tamanho da janela, desativar os nomes dos eixos, os ticks do gráfico, a barra do colourmap, etc.

printGnuplot() -> Tem como objetivo fazer o gráfico própriamente dito, traduzindo um array unidimensional em um grid 2d, produzindo primeiro o colourmap e depois o vector plot em cima

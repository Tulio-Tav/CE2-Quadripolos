# Código para Disciplina de Circuitos Elétricos II
#
# Este código simula uma linha de transmissão com três transformadores e três cargas.
# O código calcula a matriz de transmissão da linha e resolve o sistema de equações
# para encontrar a tensão e corrente em cada carga.
# O código também simula a linha com ajuste de tap e com reatores em paralelo com as cargas.
#
# Data da ultima modificação: 19-04-2025
#
#   Autores:
#   Elmer Pimentel Farias – 122110283
#   Lucas de Souza Salvino – 122210041
#   Eduardo Henrique de Freitas Coura – 122210121
#   João Henrique Morais de Nascimento – 122110054
#   Túlio Rafael de Aguiar Tavares – 120110898
#

#           +------+    Rf       Xf      +--------+                  +----------+           +--------+       
# Vac {~}---|  ~   |---/\/\/\---^^^^^----|   T1   |==LT1===||==LT3===|    T2    |===LT4=====|   T3   |-----+
#           | 60Hz |                     |69/500kV|==LT2===||-+      |500/230kV |-+         |230/69kV|     |   
#           +------+                     +--------+           |      +----------+ |         +--------+     |   
#                                                             |                   |                       [Z3]
#                                                             |                   |                        |
#                                                             |                   |                        |
#                                                            [Z1]                [Z2]                     GND     
#                                                             |                   |
#                                                            GND                 GND


import numpy as np

f  = 60                     # frequencia em Hz

w  = 2*f*np.pi              # frequencia angular em rad/s

Rf, Xf = 2, 0.38j           # parametros de impedância em serie com o transformador T1

Zth = Rf + Xf               # impedancia serie Thevenin


#-=-=-=-=-=-=-=-=-=Parametros do transformador-=-=-=-=-=-=-=-=-=#
#
# Aqui, consideramos o seguinte circuito equivalente:
#                                                        ( N1 : N2 ) 
#      ─────R1─────jX1────┬───────────R2─────jX2────────────┬  ┬─────────┐
#                         │                                 │  |        
#                       [ Rm  ]                             3||Ɛ        
#                       [ jXm ]                             3||Ɛ       
#                         │                                 3||Ɛ        
#      ───────────────────┴─────────────────────────────────┴  ┴─────────┘


R1, X1, = 7.6e-3, 3.8e-3j   # parametros gerais do transformador

R2, X2 = 33.9e-3, 0.85e-3j  # parametros gerais do transformador

Z1 = R1 + X1                # impedancia 1 do transformador

Z2 = R2 + X2                # impedancia 2 do transformador

#   impedancias "shunt" dos transformadores
#   ZT = Rm + jXm

#   Rm1 = 4320
#   Xm1 = 5050
ZT1 =  (4320 * 5050j) / (4320 + 5050j)

#   Rm2 = 432000
#   Xm2 = 505000
ZT2 = (432000 * 505000j) / (432000 + 505000j)

#   Rm3 = 402000
#   Xm3 = 607000
ZT3 = (402000 * 607000j) / (402000 + 607000j)

#Cargas

# R1 = 800ohm
# L1 = 4.1H
Zc1 = 800 + (1j* w * 4.1)

# R2 = 135.55ohm
# L2 = 0.83H
Zc2 = 135.55 + (1j* w * 0.83)

# R3 = 64.9ohm
# L3 = 0.32H
Zc3 = 64.9 + (1j * w * 0.32)

#-=-=-=-=-=-=-=-=-=Modelos de Transmissao-=-=-=-=-=-=-=-=-=#

def TransformadorIdeal(N1,N2):
    #esta funcao retorna a matriz de transmissao de um transformador ideal
    #considerando que I2 sai do pelo a direita

    matriz_Transformador = np.array([[ N1/N2 ,   0   ],
                       [   0   , N2/N1 ]])

    return matriz_Transformador

def ImpedanciaSerie(Z):

    # Matriz de impedancia

    matriz_Z = np.array([[ 1 , Z ],
                         [ 0 , 1 ]])

    return matriz_Z

def AdmitanciaShunt(Y):

    # Matriz de admitancia / Matriz inversa da impedancia

    matriz_Y = np.array([[ 1 , 0 ],
                         [ Y , 1 ]])
    return matriz_Y

def CircuitoT(Z1, Z2, Y):

    matriz_T = np.array([[ 1 + (Y * Z1) , Z1 + Z2 + (Y * Z1 * Z2) ],
                         [      Y       ,       1 + (Y * Z2)      ]])

    return matriz_T

def CircuitoPI(Z, Y1, Y2):

    matriz_Pi = np.array([[ 1 + (Y2 * Z)          ,       Z     ],
                         [ Y1 + Y2 + (Y1* Y2 * Z) , 1 + (Y1 * Z)]])

    return matriz_Pi


def LinhaDeTransmissao(Comprimento):

    #a entrada de comprimento deve ser feita em Km

    R = 0.172 * Comprimento

    L = 1j*w * 2.18e-3  * Comprimento

    C1 =  0.0136e-6 * Comprimento

    C = C1/2

    C_fasorial = 1 / (1j*120* np.pi * C)

    Z = R + L

    Y1 = 1 / C_fasorial

    Y2 = 1 / C_fasorial
    return CircuitoPI(Z, Y1, Y2)

def Cascata(*matrizes):
    resultado = matrizes[0]
    for matriz in matrizes[1:]:
        resultado = np.dot(resultado, matriz)
    return resultado


def QuadripoloParalelo(matriz1, matriz2):

    Aa, Ba, Ca, Da = matriz1[0][0], matriz1[0][1], matriz1[1][0], matriz1[1][1]

    Ab, Bb, Cb, Db = matriz2[0][0], matriz2[0][1], matriz2[1][0], matriz2[1][1]

    den = Ba + Bb

    A = (( Aa * Bb ) + (Ab * Ba ) ) / den

    B = ( Ba * Bb ) / den

    C = ( Ca + Cb + ( ( Aa - Ab )*( Db - Da ) / den ) )

    D = ( (Bb * Da ) + ( Ba * Db )) / den

    matriz_Paralelo = np.array( [[ A , B],
                                 [ C , D]])

    return matriz_Paralelo


#-=-=-=-=-=-=-=-=-=Quadripolos-=-=-=-=-=-=-=-=-=#

#impedancia em serie com a fonte

serie = ImpedanciaSerie(Zth)

#transformadores

T1 = np.dot(CircuitoT(Z1, Z2, 1/ZT1), TransformadorIdeal(69, 500))

T2 = np.dot(CircuitoT(Z1, Z2, 1/ZT2), TransformadorIdeal(500, 230))

T3 = np.dot(CircuitoT(Z1, Z2, 1/ZT3), TransformadorIdeal(230, 69))

#Linhas de transmissao

LT1 = LinhaDeTransmissao(80)

LT2 = LinhaDeTransmissao(80)

LT3 = LinhaDeTransmissao(120)

LT4 = LinhaDeTransmissao(100)

#Cargas

CargaZ1 = AdmitanciaShunt(1/Zc1)

CargaZ2 = AdmitanciaShunt(1/Zc2)

CargaZ3 = AdmitanciaShunt(1/Zc3)

#-=-=-=-=-=-=-=-=-=Simulacao da linha sem alterações-=-=-=-=-=-=-=-=-=#

ABCD = Cascata(serie, T1, QuadripoloParalelo(LT1, LT2,),
                         CargaZ1, LT3, T2, CargaZ2, LT4, T3, CargaZ3)
print('-=-'*10, 'Linha de transmissão original', '-=-'*10, '\n')

print('Matriz da linha de transmissão: \n', ABCD,'\n')

#-=-=-=-=-=-=-=-=-=solucao do sitema para Z3 -=-=-=-=-=-=-=-=-=#

A = ABCD[0][0]

B = ABCD[0][1]

C = ABCD[1][0]

D = ABCD[1][1]

Eqs = np.array([[A + (B /Zc3), 0], [-(C+(D/Zc3)), 1]])

Igualdade = np.array([69e3*np.sqrt(2), 0])

solucao = np.linalg.solve(Eqs, Igualdade)

print(f'Para a carga Z3, V = {round(np.abs(solucao[0]),4)} ∠ {np.angle(solucao[0])} V\n')

print(f'Para a carga Z3, I = {np.abs(solucao[0]/Zc3)} ∠  {np.angle(solucao[0]/Zc3)} A \n')

#-=-=-=-=-=-=-=-=-=solucao do sitema para Z2 -=-=-=-=-=-=-=-=-=#

ABCD_Z2 = Cascata(serie, T1, QuadripoloParalelo(LT1, LT2,),
                         CargaZ1, LT3, T2, CargaZ2)

A2 = ABCD_Z2[0][0]

B2 = ABCD_Z2[0][1]

C2 = ABCD_Z2[1][0]

D2 = ABCD_Z2[1][1]

Eqs2 = np.array([[A2 + (B2 /Zc2), 0], [-(C2+(D2/Zc2)), 1]])

solucao2 = np.linalg.solve(Eqs2, Igualdade)

print(f'Para a carga Z2, V = {np.abs(solucao2[0])} ∠ {np.angle(solucao2[0])} V \n')

print(f'Para a carga Z2, I = {np.abs(solucao2[0]/Zc2)} ∠  {np.angle(solucao2[0]/Zc2)} A \n') #


#-=-=-=-=-=-=-=-=-=solucao do sitema para Z1 -=-=-=-=-=-=-=-=-=#

ABCD_Z3 = Cascata(serie, T1, QuadripoloParalelo(LT1, LT2,),
                         CargaZ1)

A3 = ABCD_Z3[0][0]

B3 = ABCD_Z3[0][1]

C3 = ABCD_Z3[1][0]

D3 = ABCD_Z3[1][1]

Eqs3 = np.array([[A3 + (B3 /Zc1), 0], [-(C3+(D3/Zc1)), 1]])

solucao3 = np.linalg.solve(Eqs3, Igualdade)

print(f'Para a carga Z1, V = {np.abs(solucao3[0])} ∠ {np.angle(solucao3[0])} V \n')

print(f'Para a carga Z1, I = {np.abs(solucao3[0]/Zc1)} ∠  {np.angle(solucao3[0]/Zc1)} A \n')


#-=-=-=-=-=-=-=-=-=Simulacao com ajuste de tap -=-=-=-=-=-=-=-=-=#



T1_tap = np.dot(CircuitoT(Z1, Z2, 1/ZT1), TransformadorIdeal(69, 350.1092))

T2_tap = np.dot(CircuitoT(Z1, Z2, 1/ZT2), TransformadorIdeal(500, 228.664))

T3_tap = np.dot(CircuitoT(Z1, Z2, 1/ZT3), TransformadorIdeal(230, 67.2282))

ABCD_tap = Cascata(serie, T1_tap, QuadripoloParalelo(LT1, LT2,),
                         CargaZ1, LT3, T2_tap, CargaZ2, LT4, T3_tap, CargaZ3)

print('-=-'*10, 'Linha de transmissão com ajuste de tap', '-=-'*10, '\n')

print('Matriz da linha de transmissão: \n', ABCD_tap,'\n')

A_tap = ABCD_tap[0][0]

B_tap = ABCD_tap[0][1]

C_tap = ABCD_tap[1][0]

D_tap = ABCD_tap[1][1]

Eqs_tap = np.array([[A_tap + (B_tap /Zc3), 0], [-(C_tap+(D_tap/Zc3)), 1]])

solucao_tap = np.linalg.solve(Eqs_tap, Igualdade)

print(f'Para a carga Z3 com ajuste de tap, V = {np.abs(solucao_tap[0])} ∠ {np.angle(solucao_tap[0])} V \n')

print(f'Para a carga Z3 com ajuste de tap, I = {np.abs(solucao_tap[0]/Zc3)} ∠  {np.angle(solucao_tap[0]/Zc3)} A \n')

#-=-=-=-=-=-=-=-= Z2 tap -=-=-=-=-=-=-=-=-=-=-=-=

ABCD_tap_Z2 = Cascata(serie, T1_tap, QuadripoloParalelo(LT1, LT2,),
                         CargaZ1, LT3, T2_tap, CargaZ2)

A_tap_Z2 = ABCD_tap_Z2[0][0]

B_tap_Z2 = ABCD_tap_Z2[0][1]

C_tap_Z2 = ABCD_tap_Z2[1][0]

D_tap_Z2 = ABCD_tap_Z2[1][1]

Eqs_tap_Z2 = np.array([[A_tap_Z2 + (B_tap_Z2 /Zc2), 0], [-(C_tap_Z2+(D_tap_Z2/Zc2)), 1]])

solucao_tap_Z2 = np.linalg.solve(Eqs_tap_Z2, Igualdade)

print(f'Para a carga Z2 com ajuste de tap, V = {np.abs(solucao_tap_Z2[0])} ∠ {np.angle(solucao_tap_Z2[0])} V \n')

print(f'Para a carga Z2 com ajuste de tap, I = {np.abs(solucao_tap_Z2[0]/Zc2)} ∠  {np.angle(solucao_tap_Z2[0]/Zc2)} A \n')

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-= Z1 tap -=-=-=-=-=-=-=-=-=

ABCD_tap_Z1 = Cascata(serie, T1_tap, QuadripoloParalelo(LT1, LT2,),
                         CargaZ1)

A_tap_Z1 = ABCD_tap_Z1[0][0]

B_tap_Z1 = ABCD_tap_Z1[0][1]

C_tap_Z1 = ABCD_tap_Z1[1][0]

D_tap_Z1 = ABCD_tap_Z1[1][1]

Eqs_tap_Z1 = np.array([[A_tap_Z1 + (B_tap_Z1 /Zc1), 0], [-(C_tap_Z1+(D_tap_Z1/Zc1)), 1]])

solucao_tap_Z1 = np.linalg.solve(Eqs_tap_Z1, Igualdade)

print(f'Para a carga Z1 com ajuste de tap, V = {np.abs(solucao_tap_Z1[0])} ∠ {np.angle(solucao_tap_Z1[0])} V \n')

print(f'Para a carga Z1 com ajuste de tap, I = {np.abs(solucao_tap_Z1[0]/Zc1)} ∠  {np.angle(solucao_tap_Z1[0]/Zc1)} A \n')

#-=-=-=-=-=-=-=-=-=Simulacao com reatores -=-=-=-=-=-=-=-=-=#

reator1 = AdmitanciaShunt( 1/(1j * w * 650e-3))

reator2 = AdmitanciaShunt( 1/(0.001 + 1j * w * 3))

reator3 = AdmitanciaShunt( 1/(0.001+ 1j * w * 1.25))

ABCD_reator = Cascata(serie, T1, QuadripoloParalelo(LT1, LT2,),
                         CargaZ1, reator1, LT3, T2, CargaZ2, reator2, LT4, T3, CargaZ3, reator3)

print('-=-'*10, 'Linha de transmissão com associação de reatores em paralelo com a carga', '-=-'*10, '\n')

print('Matriz da linha de transmissão: \n', ABCD_reator,'\n')

A_reator = ABCD_reator[0][0]

B_reator = ABCD_reator[0][1]

C_reator = ABCD_reator[1][0]

D_reator = ABCD_reator[1][1]

Eqs_reator = np.array([[A_reator + (B_reator /Zc3), 0], [-(C_reator+(D_reator/Zc3)), 1]])

Igualdade_reator = np.array([69e3*np.sqrt(2), 0])

solucao_reator = np.linalg.solve(Eqs_reator, Igualdade_reator)

print(f'Para a carga Z3 com reator em paralelo, V = {round(np.abs(solucao_reator[0]),4)} ∠ {np.angle(solucao_reator[0])} V \n')

print(f'Para a carga Z3 com reator paralelo, I = {np.abs(solucao_reator[0]/Zc3)} ∠  {np.angle(solucao_reator[0]/Zc3)} A \n')

ABCD_reator_Z2 = Cascata(serie, T1, QuadripoloParalelo(LT1, LT2,),
                         CargaZ1, reator1, LT3, T2, CargaZ2, reator2)

A_reator_Z2 = ABCD_reator_Z2[0][0]

B_reator_Z2 = ABCD_reator_Z2[0][1]

C_reator_Z2 = ABCD_reator_Z2[1][0]

D_reator_Z2 = ABCD_reator_Z2[1][1]

Eqs_reator_Z2 = np.array([[A_reator_Z2 + (B_reator_Z2 /Zc2), 0], [-(C_reator_Z2+(D_reator_Z2/Zc2)), 1]])

Igualdade_reator_Z2 = np.array([69e3*np.sqrt(2), 0])

solucao_reator_Z2 = np.linalg.solve(Eqs_reator_Z2, Igualdade_reator_Z2)

print(f'Para a carga Z2 com reator em paralelo, V = {round(np.abs(solucao_reator_Z2[0]),4)} ∠ {np.angle(solucao_reator_Z2[0])} V \n')

print(f'Para a carga Z2 com reator paralelo, I = {np.abs(solucao_reator_Z2[0]/Zc2)} ∠  {np.angle(solucao_reator_Z2[0]/Zc2)} A \n')

ABCD_reator_Z1 = Cascata(serie, T1, QuadripoloParalelo(LT1, LT2,),
                         CargaZ1, reator1)

A_reator_Z1 = ABCD_reator_Z1[0][0]

B_reator_Z1 = ABCD_reator_Z1[0][1]

C_reator_Z1 = ABCD_reator_Z1[1][0]

D_reator_Z1 = ABCD_reator_Z1[1][1]

Eqs_reator_Z1 = np.array([[A_reator_Z1 + (B_reator_Z1 /Zc2), 0], [-(C_reator_Z1+(D_reator_Z1/Zc1)), 1]])

Igualdade_reator_Z1 = np.array([69e3*np.sqrt(2), 0])

solucao_reator_Z1 = np.linalg.solve(Eqs_reator_Z1, Igualdade_reator_Z1)

print(f'Para a carga Z1 com reator em paralelo, V = {round(np.abs(solucao_reator_Z1[0]),4)} ∠ {np.angle(solucao_reator_Z1[0])} V \n')

print(f'Para a carga Z1 com reator paralelo, I = {np.abs(solucao_reator_Z1[0]/Zc1)} ∠  {np.angle(solucao_reator_Z1[0]/Zc1)} A \n')





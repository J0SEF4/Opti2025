from gurobipy import Model, GRB, quicksum
import numpy as np
import pandas as pd

# datos/parámetros del modelo

#   conjuntos
T= #meses de planificación, t
R= #relaves, r
F= #fuentes de agua, f
N= #nodos de la red, fuentes, nodos intermedios y relaves, n
A= #arcos de la red, (i,j)

#   parámetros
C_agua = #costo por ton de agua aplicado en relave r, mes t
C_veginst =  #costo instalar cubierta vegetal en relave r, mes t
C_vegmant = #costo mantención cubierta vegetal en relave r, mes t
C_mant = #costo oportunidad realizar mantención en relave r, mes t
PM_base = #concentracion inicial de PM en relave r
PM_agregado = #concentracion de PM que se agrega al relave r, mes t
PM_max = #umbral maximo de PM permitido en relave r por mes de planificación t
PM = #umbral promedio maximo de PM permitido para todos los relaves r en mes t
alpha = #reducción PM por cada ton de agua aplicada en relave r
beta = #maxima reducción PM por cubierta vegetal en relave r por mes de planificación
H_mant = #agua minima necesaria para humedecer relave r cuando se realiza mantención por mes de planificación
H_cubierta = #agua minima necesaria para humedecer una cubierta vegetal a su nivel optimo en relave r por mes de planificación
B_max = #presupuesto maximo a lo largo del periodo de planificación
P = #periodos de duración de permanencia de cubierta vegetal en relave r 
P_max = #cantidad maxima de periodos consecutivos que pueden pasar sin hacer mantención a un relave
P_cubierta = #cantidad  maxima de periodos consecutivos que pueden pasar sin satisfacer completamente el requerimiento de agua de la cubierta
W_base = #cantidad agua inicial en fuente f
C_flujo = #costo transportar 1 ton de agua por arco (i.j)
K_flujo = #capacidad maxima de flujo de agua en argo (i,j) mensualmente
U_agua = #cantidad agua maxima que se puede aplicar al relave r en cualquier periodo sin considerar cubierta vegetal
a = #ton de agua minima que se aplica a un relave r si se decide humedecerlo
W_entrante = #flujo neto entrante a fuente f durante mes t.

# Generar el modelo
model = Model()

#variables
x_agua = model.addVars(R, T, vtype=GRB.CONTINUOUS, name="x_agua", lb=0)              
x_cubierta = model.addVars(R, T, vtype=GRB.CONTINUOUS, name="x_cubierta", lb=0)      
x_flujo = model.addVars(A, T, vtype=GRB.CONTINUOUS, name="x_flujo", lb=0)           

y_veg = model.addVars(R, T, vtype=GRB.BINARY, name="y_veg")
y_mant = model.addVars(R, T, vtype=GRB.BINARY, name="y_mant")
y_agua = model.addVars(R, T, vtype=GRB.BINARY, name="y_agua")

z_veg = model.addVars(R, T, vtype=GRB.BINARY, name="z_veg")

r = model.addVars(R, T, vtype=GRB.BINARY, name="r")

PM = model.addVars(R, T, vtype=GRB.CONTINUOUS, name="PM", lb=0)               
W = model.addVars(F, T, vtype=GRB.CONTINUOUS, name="W", lb=0)                 
beta_efectivo = model.addVars(R, T, vtype=GRB.CONTINUOUS, name="beta_efectivo", lb=0)        

# Llamamos a update
model.update()

#restricciones
#R1 Consistencia de agua
for r in R:
    for t in T:
        model.addConstr(a[r]*y_agua[r,t] <= x_agua[r,t])
        model.addConstr(x_agua[r,t] <= U_agua[r]*y_agua[r,t])

#R2 Capacidad de cubierta
for r in R:
    for t in T:
        model.addConstr(x_cubierta[r,t] <= H_cubierta[r]*z_veg[r,t])

#R3 Restistencia de cubierta - REVISAR
for r in R:
    for t in T:
        model.addConstr(quicksum(1-r[r,t2] for t2 in range(t,min(t+P_cubierta,T)))) <= (P_cubierta*z_veg[r,t])+(T*(1 - z_veg[r,t]))

#R4 Reducción de PM efectiva
for r in R:
    for t in T:
        model.addConstr(beta_efectivo[r,t] <= (x_cubierta[r,t]/H_cubierta[r])*beta[r])
        model.addConstr(beta_efectivo[r,t] <= z_veg[r,t]*beta[r])

#R5 Flujo de la concentración de material particulado
for r in R:
    #mes=1
    model.addConstr(PM_base[r] + PM_agregado[r,0] - alpha[r]*(x_agua[r,0] - x_cubierta[r,0]) - beta_efectivo[r,0] == PM[r,0])
    #meses posteriores
    for t in T:
        if t > 0: #como empezamos en t=0, t=1 es el 2do mes
            model.addConstr(PM[r,t-1] + PM_agregado[r,t] - alpha[r]*(x_agua[r,t] - x_cubierta[r,t]) - beta_efectivo[r,t] == PM[r,t])

#R6 Control de PM por relave
for r in R:
    for t in T:
        model.addConstr(PM[r,t] <= PM_max[r,t])

#R7 Control de PM promedio
for t in T:
    model.addConstr((1/len(R))*quicksum(PM[r,t] for r in R) <= PM[t])

#R8 Agua mínima si hay mantención
for r in R:
    for t in T:
        model.addConstr(x_agua[r,t] >= H_mant[r]*y_mant[r,t])

#R9 Ventana de mantención obligatoria
for r in R:
    for t in T:
        model.addConstr(quicksum(y_mant[r, t2] for t2 in range(max(0, t-P_max+1), t+1)) >= 1)

#R10 Acciones con mantención activa
for r in R:
    for t in T:
        model.addConstr(y_veg[r,t] <= y_mant[r,t])

#R11 Presencia de la cubierta vegetal

#R12 Presupuesto total
model.addConstr(
    quicksum(C_agua[r,t]*x_agua[r,t] + C_vegmant[r,t]*z_veg[r,t] + C_mant[r,t]*y_mant[r,t] + C_veginst[r,t]*y_veg[r,t] for r in R for t in T) 
    + quicksum(C_flujo[(i,j)]*x_flujo[(i,j),t] for (i,j) in A for t in T)
    <= B_max
)
#R13 Balance de flujo en nodos intermedios

#R14 Condición inicial de la fuente
for f in F:
    model.addConstr(W[f,0] == W_base[f])

#R15 Balance de la fuente para cada mes t 

#R16 No extraer más de lo disponible

#R17 Nodos de relave
for r in R:
    for t in T:
        model.addConstr(quicksum(x_flujo[(j,r),t] for (j,r) in A) == x_agua[r,t])
        
#R18 Capacidad de las tuberías
for (i, j) in A:
    for t in T:
        model.addConstr(x_flujo[(i,j), t] <= K_flujo[A.index((i, j))])


#función objetivo
model.setObjective(quicksum(PM[r,t] for r in R for t in T), GRB.MINIMIZE)

#resolver modelo
model.optimize()

#imprimir resultados
if model.Status == GRB.OPTIMAL:
    print("\nValor óptimo de la función objetivo:", model.ObjVal)
    for r in R:
        for t in T:
            print(f"Relave {r}, Mes {t}, PM = {PM[r, t].X:.2f}, xagua = {x_agua[r, t].X:.2f}")
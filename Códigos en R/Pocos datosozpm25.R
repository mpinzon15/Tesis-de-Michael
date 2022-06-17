library(R2jags)
library(coda)
library(boa)
library(lattice)
rm(list=ls(all=TRUE))
excedenciaoz <-read.csv2("O3_generada_240421.csv", header=T, sep = ",")
excedenciapm10 <-read.csv2("PM10_generada_240421.csv", header=T, sep = ",")
excedenciapm25 <-read.csv2("PM25_generada_240421.csv", header=T, sep = ",") 
v_exc=NULL
v_exc1=NULL
v_exc2=NULL
for (i in 1:500){if (excedenciaoz$x[i]>=50.9375){v_exc2<-c(v_exc2,i)}}
for (i in 1:500){if (excedenciapm10$x[i]>=75){v_exc1<-c(v_exc1,i)}}
for (i in 1:500){if (excedenciapm25$x[i]>=37){v_exc<-c(v_exc,i)}}

a<-v_exc2[1:100] #Ozono
b<-v_exc1[1:100] #pm10
c<-v_exc[1:100]  #pm2.5


#exced_03<-read.csv("excedencias_oz.csv")
#exced_pm10<-read.csv("excedencias_pm10.csv")
#exced_pm25<-read.csv("excedencias_pm25.csv")

d1<-v_exc2
d1_0<-c(0,d1)
n=length(d1)
T=500


d2<-v_exc
d2_0<-c(0,d2)
K=length(d2) 

dif_pm10=rep(0,length(d1))
for(i in 1:length(dif_pm10)){
  dif_pm10[i]=d1_0[i+1]-d1_0[i]
}
d1_dif=dif_pm10


dif_pm25=rep(0,length(d2))
for(i in 1:length(dif_pm25)){
  dif_pm25[i]=d2_0[i+1]-d2_0[i]
}
d2_dif=dif_pm25


mod_string = "
data {
for(i in 1:n){
for(j in 1:K){
dummy[i,j]<- 0
}
}}
model{
for(j in 1:K){
for(i in 1:n){
phi3[i,j] <- -log(L[i,j])
dummy[i,j] ~ dpois(phi3[i,j])
log(lambda[i,j])<- log((alpha[1]/sigma[1])*pow(d1_dif[i]/sigma[1],alpha[1]-1)
*(alpha[2]/sigma[2])*pow(d2_dif[j]/sigma[2],alpha[2]-1)
*(1+theta*(1-2*exp(-pow(d2_dif[j]/sigma[2],alpha[2])))
*(1-2*exp(-pow(d1_dif[i]/sigma[1],alpha[1]))))
/(1+theta*(1-exp(-pow(d1_dif[i]/sigma[1],alpha[1])))
*(1-exp(-pow(d2_dif[j]/sigma[2],alpha[2])))))
L[i,j]<-lambda[i,j]*F
m[i,j]<-pow(d1[i]/sigma[1],alpha[1])
+pow(d2[j]/sigma[2],alpha[2])-log(1
+theta*(1-exp(-pow(d1[i]/sigma[1],alpha[1])))
*(1-exp(-pow(d2[j]/sigma[2],alpha[2]))))
lambda2[i,j]<- (alpha[1]/sigma[1])*pow(d1[i]/sigma[1],alpha[1]-1)
*(alpha[2]/sigma[2])*pow(d2[j]/sigma[2],alpha[2]-1)
*(1+theta*(1-2*exp(-pow(d2[j]/sigma[2],alpha[2])))
*(1-2*exp(-pow(d1[i]/sigma[1],alpha[1]))))
/(1+theta*(1-exp(-pow(d1[i]/sigma[1],alpha[1])))
*(1-exp(-pow(d2[j]/sigma[2],alpha[2]))))
}
}
F<-exp(-(pow(T/sigma[1],alpha[1])*pow(T/sigma[2],alpha[2]))
/(n*K))
alpha[1] ~ dgamma(1627,1344)
sigma[1] ~ dgamma(9,7)
alpha[2] ~ dgamma(391,313)
sigma[2] ~ dgamma(171,16)
theta ~ dunif(-1,1)
} 
"
#alpha[1] ~ dgamma(1,4)
#sigma[1] ~ dgamma(1,12)
#alpha[2] ~ dgamma(1,4)
#sigma[2] ~ dgamma(1,12)



bivariado_oz_pm25 =list(n=n, K=K, T=T, d1_dif=d1_dif, d2_dif=d2_dif, d1=d1, d2=d2)
bivariado_oz_pm25params= c("theta","sigma","alpha","m","lambda")
modbiv_oz_pm25 = jags.model(textConnection(mod_string), data=bivariado_oz_pm25, n.chains=2)
oz_pm25.upd<-update(modbiv_oz_pm25,1000)
mod_oz_pm25_sim =coda.samples(model=modbiv_oz_pm25, bivariado_oz_pm25params, n.iter=5000, thin=10)
summary(mod_oz_pm25_sim[c(1,2)][,c(1:2,31683:31685)])$statistics[1:5,1]#500
summary(mod_oz_pm25_sim[c(1,2)][,c(1:2,12015:12017)])$statistics[1:5,1]#400
#summary(mod_oz_pm25_sim[c(1,2)][,c(1:2,20003:20005)])$statistics[1:5,1]

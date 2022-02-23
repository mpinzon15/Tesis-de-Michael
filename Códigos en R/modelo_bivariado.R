library(R2jags)
library(coda)
library(boa)
library(lattice)

rm(list=ls(all=TRUE))
ex_03<-c( 12,21,25,28,29,34,35,45,46,73,75,88,89,90,106,269,315,352,353,354,355,381,382,392,414,415,416,417,418,420,447,448,449,487,489)
excpm10<-c(10,  11,  12,  13,  14,  15,  16,  17,
           21,  22,  23,  24  ,27 , 28,  29,  30,
           31,  32,  33,  34  ,35  ,36  ,37,  38,
           39,  41,  42,  44,  45,  47,  48,  50,
           51,  52,  54,  55,  56,  57,  58,  64,
           65,  66,  70,  71,  72,  77,  78,  79,
           80,  81,  86,  87,  88,  89,  90,  91,
           92,  99, 100, 101, 102, 104, 111, 113,
           115, 116, 120, 122, 123, 125, 126, 127,
           132, 133, 134, 135, 136, 137, 139, 140,
           141, 143, 146, 147, 148, 149, 150, 154,
           155, 156, 160, 162, 169, 177, 178, 179,
           184, 190, 210, 216, 217, 223, 227, 232,
           233, 234, 240, 241, 242, 245, 248, 260,
           272, 273, 281, 283, 284, 287, 288, 289,
           293, 294, 296, 297, 301, 302, 303, 304,
           305, 308, 309, 310, 311, 312, 317, 318,
           319, 322, 323, 324, 325, 326, 328, 329,
           330, 331, 332, 333, 335, 336, 337, 338,
           339, 340, 341, 342, 343, 344, 347, 349,
           350, 351, 352, 353, 354, 356, 357, 358,
           361, 379, 385, 388, 394, 395, 398, 400,
           402, 403, 406, 407, 408, 410, 412, 413,
           414, 415, 416, 417, 419, 420, 421, 422,
           423, 424, 428, 429, 430, 436, 441, 448,
           449, 450, 451, 452, 454, 455, 456, 461,
           462, 463, 464, 465, 476, 477, 478, 479,
           485, 486, 487, 489, 490, 491, 528, 532,
           546, 547, 548, 549, 553, 555, 577, 596,
           598, 603, 604)

d1<-ex_03
d1_0<-c(0,d1)
d2<-excpm10
d2_0<-c(0,d1)
d1_dif<-(c(1:length(d1)))
d2_dif<-(c(1:length(d2)))
n=length(d1)
K=length(d2) 
T=608

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
alpha[1] ~ dunif(1,3)
sigma[1] ~ dunif(11,20)
alpha[2] ~ dunif(0.5,1.5)
sigma[2] ~ dunif(0,0.5)
theta ~ dunif(-1,1)
} 
"
bivariado_oz_pm10 =list(n=n, K=K, T=T, d1_dif=d1_dif, d2_dif=d2_dif, d1=d1, d2=d2)
bivariado_oz_pm10params= c("theta","sigma","alpha")
modbiv_oz_pm10 = jags.model(textConnection(mod_string), data=bivariado_oz_pm10, n.chains=2)
oz_pm10.upd<-update(modbiv_oz_pm10,1000)
mod_oz_pm10_sim =coda.samples(model=modbiv_oz_pm10, bivariado_oz_pm10params, n.iter=5000, n.thin=10)
summary(mod_oz_pm10_sim[c(1,2)][,c(1:2,497:500)])$statistics[1:5,1]
#alpha[1]    alpha[2]    sigma[1]     sigma[2]   theta 
#1.00004342  0.51777298  16.62158737  0.01881131 -0.99992915
#alpha[1]    alpha[2]    sigma[1]    sigma[2]       theta 
#1.00004376  0.51783197 14.94230262  0.02342856 -0.99994252 
summary(mod_oz_pm10_sim)
summary(mod_oz_pm10_sim[c(1,2)][,c(1:2,6997:7000)])$statistics[1:5,2]
#alpha[1]     alpha[2]     sigma[1]     sigma[2]        theta 
#4.251542e-05 5.745281e-03 2.385305e+00 6.648830e-03 7.492540e-05

library(sp)
rm(list=ls(all=TRUE))
library(readr)
excedenciapm25 <-read.csv2("PM25_generada_240421.csv", header=T, sep = ",") 
i<- NULL; v_exc<-NULL
excedenciapm25$x <-as.numeric(as.character(excedenciapm25$x))
#mean(excedenciapm25$x)
#var(excedenciapm25$x)
for (i in 1:1096){if (excedenciapm25$x[i]>=37){v_exc<-c(v_exc,i)}}

pm25f<- cbind(mis_fechas,excedenciapm25)#NUEVO
pm25fil<-pm25f %>% filter(x>=37)


library(R2jags)
library(coda)
library(boa)
library(lattice)
d=v_exc
K=length(d)
R=length(excedenciapm25$x)
mod_string = "
data{for(i in 1:K){
dummy[i] <- 0}
}
model{
for (i in 1:K) {
phi[i] <- -log(L[i])
dummy[i] ~ dpois(phi[i])
log(lambda[i]) <-  log( (alpha/sigma) * pow(d[i]/sigma,alpha-1) )
L[i] <- lambda[i]*F 
m[i]<-pow(d[i]/sigma,alpha)
}
F <-  exp(-pow(R/sigma,alpha)/K)
alpha ~dunif(1,2.8)
sigma ~ dunif(10, 100) 
} 
"

#alpha ~ dunif(1,3.5)
#sigma ~ dunif(10, 80)
pm2.5_sinpc_45.data =list(d=d, K=K, R=R)
pm2.5_sinpc_45.params= c("alpha","sigma","m")

mod_pm2.5_sinpc_45= jags.model(textConnection(mod_string), data=pm2.5_sinpc_45.data, n.chains=5)
pm2.5_sinpc_45.upd<-update(mod_pm2.5_sinpc_45, 11000)# burn-in
#dic0  <- dic.samples(mod_pm2.5_sinpc_45, variable.names=pm2.5_sinpc_45.params, n.iter=20000, progress.bar="none")
mod_pm2.5_sinpc_45_sim =coda.samples(model=mod_pm2.5_sinpc_45,variable.names=pm2.5_sinpc_45.params,n.iter=20000,n.thin=15)
#pr_pm2.5_sinpc_45<-print(mod_pm2.5_unpc_45_sim)

m_sim0<-summary(mod_pm2.5_sinpc_45_sim[c(1:5)][,c(2:333)])$statistics[1:332,1]
m2.5p_0<-summary(mod_pm2.5_sinpc_45_sim[c(1,2)][,c(2:333)])$quantiles[1:332,1]
m97.5p_0<-summary(mod_pm2.5_sinpc_45_sim[c(1,2)][,c(2:333)])$quantiles[1:332,5]


pdf("traza1a.pdf")

plot(mod_pm2.5_sinpc_45_sim[c(1:5)][,c(1,334)])
Xdev.off()
densityplot(mod_pm2.5_sinpc_45_sim[c(1:5)][,c(1,334)])
pdf("densidad1.pdf")
xyplot(mod_pm2.5_sinpc_45_sim[c(1:5)][,c(1,334)])
dev.off()

autocorr.plot(mod_pm2.5_sinpc_45_sim[c(1:5)][,c(1,334)],100)
gelman.plot(mod_pm2.5_sinpc_45_sim[c(1:5)][,c(1,334)])
gelman.diag(mod_pm2.5_sinpc_45_sim[c(1:5)][,c(1,334)])
geweke.plot(mod_pm2.5_sinpc_45_sim[c(1:5)][,c(1,334)])
geweke.diag(mod_pm2.5_sinpc_45_sim[c(1:5)][,c(1,334)])
raftery.diag(mod_pm2.5_sinpc_45_sim[c(1:5)][,c(1,334)])
heidel.diag(mod_pm2.5_sinpc_45_sim[c(1:5)][,c(1,334)])



s1_sinpc_45<-summary(mod_pm2.5_sinpc_45_sim[c(1:5)][,c(1,334)])$statistics[1:2,1]
s2_sinpc_45<-summary(mod_pm2.5_sinpc_45_sim[c(1:5)][,c(1,334)])$statistics[1:2,2]
s3_sinpc_45<-summary(mod_pm2.5_sinpc_45_sim[c(1:5)][,c(1,334)])$quantiles[1:2,1]
s4_sinpc_45<-summary(mod_pm2.5_sinpc_45_sim[c(1:5)][,c(1,334)])$quantiles[1:2,5]

summary_sinpc_45 <- matrix(c(s1_sinpc_45[1], s1_sinpc_45[2], 
                             s2_sinpc_45[1], s2_sinpc_45[2], 
                             s3_sinpc_45[1], s3_sinpc_45[2],
                             s4_sinpc_45[1], s4_sinpc_45[2]), c(2, 4))
rownames(summary_sinpc_45) <- c("$\\alpha$", "$\\sigma$") 
colnames(summary_sinpc_45) <- c("mean", "sd","$2.5\\%$","$97.5\\%$") 
library("xtable")
summary_sinpc_45 <- xtable(summary_sinpc_45) 
print(summary_sinpc_45, sanitize.text.function = function(x) {x})

dic0_45  <- dic.samples(mod_pm2.5_sinpc_45, variable.names=pm2.5_sinpc_45.params, n.iter=20000, progress.bar="none")
vec_dic0_45<-c(round(sum(dic0_45$deviance),0), round(sum(dic0_45$penalty),4), round(sum(dic0_45$deviance)+sum(dic0_45$penalty)))
####Graficas##############################################################

par(mfrow=c(2,1))
par(mai=c(0.9,0.8,0.2,0.2))
plot(x = d, y = 1:length(d),yaxt='n',type = 'l',xlab = 'Days',ylab = 'means_PM2.5_spc',pch = 3)
lines(d,m_sim0,col='blue', type="l", lty=2,lwd = 2) 
lines(d, m2.5p_0, yaxt='n', type = 'l', xlab = 'Days',ylab = 'means_PM2.5_45',pch = 1, col='green', lwd = 1)
lines(d, m97.5p_0, yaxt='n', type = 'l', xlab = 'Days',ylab = 'means_PM2.5_45',pch = 1, col='green', lwd = 1)
axis(2, at = seq(0, 300, 20), las=2)

pm25fil

plot(x = pm25fil$mis_fechas, y = 1:length(pm25fil$mis_fechas),xaxt='n',type = 'l',xlab = 'Days',ylab = 'means_PM2.5',pch = 3)
axis(1, 
     pm25fil$mis_fechas, format(pm25fil$mis_fechas, "%b/%y"), las=2, cex.axis=0.72)
lines(pm25fil$mis_fechas,m_sim0,col='blue', type="l", lty=2,lwd = 2) 
lines(pm25fil$mis_fechas, m2.5p_0, yaxt='n', type = 'l', xlab = 'Days',ylab = 'means_PM2.5',pch = 1, col='green', lwd = 1)
lines(pm25fil$mis_fechas, m97.5p_0, yaxt='n', type = 'l', xlab = 'Days',ylab = 'means_PM2.5',pch = 1, col='green', lwd = 1)




 # DIC SIN PUNTOS DE CAMBIO ------------------------------------------------

DIC_sinpc <- matrix(c(vec_dic0_45[1], vec_dic0_45[2], vec_dic0_45[3],
                      vec_dic0_50[1], vec_dic0_50[2], vec_dic0_50[3]), c(2, 3),  byrow=TRUE)
rownames(DIC_sinpc) <- c("$PM2.5$ Umbral  45", "$PM2.5$ Umbral 50") 
colnames(DIC_sinpc) <- c("Mean deviance", "penalty", "Penalized deviance") 
library("xtable")
DIC_sinpc <- xtable(DIC_sinpc) 
print(DIC_sinpc, sanitize.text.function = function(x){x})
  

# GRAFICAS CON LOS PARAMETROS 

funcionpm25 <- function(x){(a2/b2)*(x/b2)^(a2-1)}
x=seq(1,1096)
a2<-1.25
b2<-10.87

plot(x=mis_fechas, y=funcionpm25(x), type="l", xaxt='n', xlab="Days", ylab="Rate_function_pm25")
axis(1, 
     mis_fechas, format(mis_fechas, "%b/%y"), las=2, cex.axis=0.72)

#plot(funcionpm25(x), type="l", xlab="Days", ylab="Rate_function_pm25")

?dgamma


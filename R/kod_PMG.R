# vari to zmienna panelowa np. vari=dane[,5]
# quan to wektor licznosci jednostek w panelu
# zaklada sie, ze jednostki nie musza byc jednakowo liczne T1<>T2<>T3...<>Tn
# postac wektora: quan=c(10,5,8), T1=10, T2=5, T3=8
DiffPanel=function(variable, quantity)
{ # pocz function
# adresy poczatkow i koncow zmiennych jednostek panelu
vari=variable
quan=quantity
do=cumsum(quan)
od=(do-quan)+1
if (is.data.frame(vari)) variDIFF=as.matrix(as.numeric(vari[[1]])) else
variDIFF=as.matrix(as.numeric(vari))
DIFF=NULL
for (i in 1:length(quan))
{ # pocz i
pomDIFF=as.numeric(variDIFF[od[i]:do[i],1])
pomDIFF=c(NA,(pomDIFF[2:length(pomDIFF)]-pomDIFF[1:(length(pomDIFF)-1)]))
DIFF=c(DIFF,pomDIFF)
} # kon i
DiffPanelPom=as.matrix(DIFF)
if (!is.null(row.names(vari))) row.names(DiffPanelPom)=row.names(vari)
DiffPanel=DiffPanelPom
} # kon function

# vari to zmienna panelowa np. vari=dane[,5]
# quan to wektor licznosci jednostek w panelu
# zaklada sie, ze jednostki nie musza byc jednakowo liczne T1<>T2<>T3...<>Tn
# postac wektora: quan=c(10,5,8), T1=10, T2=5, T3=8
LagPanel=function(variable, quantity)
{ # pocz function
vari=variable
quan=quantity
if (is.data.frame(vari)) variLAG=as.matrix(as.numeric(vari[,1])) else
variLAG=as.matrix(as.numeric(vari))
# adresy poczatkow i koncow zmiennych jednostek panelu
do=cumsum(quan)
od=(do-quan)+1
LAG=NULL
for (i in 1:length(quan))
{ # pocz i
pomLAG=as.numeric(variLAG[od[i]:do[i],1])
pomLAG=c(NA,pomLAG[1:(length(pomLAG)-1)])
LAG=c(LAG,pomLAG)
} # kon i
LagPanelPom=as.matrix(LAG)
if (!is.null(row.names(vari))) row.names(LagPanelPom)=row.names(vari)
LagPanel=LagPanelPom
} # kon function

# dataset to zbior danych panelowych w postaci stakowanych szeregow czasowych
# guan to wektor c(licznosc1, licznosc2,....,licznoscN) dlugodci szeregow czasowych
# w poszczegolnych grupach
PanelNaOmit=function(dataset, quantity)
{ # pocz function
quan=quantity
# adresy poczatkow i koncow zmiennych jednostek panelu
do=cumsum(quan)
od=(do-quan)+1
# nowy wektor licznosci T1,...,Tn dla poszczegolnych jednostek
# usuniecie obiektow NA spowoduje zmiany licznosci
nquan=NULL
# tu beda zbierane dane na.omit
P1NaOmit=NULL
for (i in 1:length(quan))
{ # pocz i
P2NaOmit=na.omit(dataset[od[i]:do[i],])
P1NaOmit=rbind(P1NaOmit,P2NaOmit)
nquan=c(nquan,nrow(P2NaOmit))
} # kon i
PanelNaOmit=list(dataset=P1NaOmit,quantity=nquan)
} # kon function

# test autokorelacji Breuscha-Godfreya
# reszty to szereg reszt oryginalnej relacji gdzie badamy autokorelacje
# Xorg to macierz zmiennych objasniajacych oryginalnej relacji
# p1 to rzad badanej autokorelacji
BGtest=function(residuals, explvariab, acor.ord)
{ # pocz test BG
reszty=residuals
Xorg=explvariab
p1=acor.ord
# zmienne objasniane relacji oryginalnej
reszty=as.numeric(reszty)
pomX_testBG=as.matrix(Xorg)
p1=as.numeric(p1)
AR_testBG=NULL
for (i in 1:p1)
{ # pocz i
pomAR=rbind(matrix(NA,i,1),as.matrix(reszty[1:(length(reszty)-i)]))
AR_testBG=cbind(AR_testBG, pomAR)
} # kon i
# wyrownanie dlugosci zmiennych (z powodu lagow)
YX_testBG=cbind(as.matrix(reszty), pomX_testBG, AR_testBG)
YX_testBG=na.omit(YX_testBG)
ModURSS=lm(YX_testBG[,1]~YX_testBG[,-1])
URSS=sum((resid(ModURSS))^2)
ModRRSS=lm(YX_testBG[,1]~YX_testBG[,2:(ncol(pomX_testBG)+1)])
RRSS=sum((resid(ModRRSS))^2)
stat_F=((RRSS-URSS)/URSS)*((nrow(YX_testBG)-ncol(YX_testBG)-1-p1)/p1)
prob_F=1-pf(stat_F, p1, ((nrow(YX_testBG)-ncol(YX_testBG)-1-p1)))
wynikF=matrix(c(stat_F, prob_F),1,2)

#stat_Chi=(T*stat_F*p1)/(T-k-1-p+stat_F*p1)
stat_Chi=(stat_F*p1)
prob_Chi=1-pchisq(stat_Chi, p1)
wynikChi=matrix(c(stat_Chi, prob_Chi),1,2)

TestBGpom=data.frame(rbind(wynikF, wynikChi), row.names=c("F stat","Chi^2 stat"))
names(TestBGpom)<-c("stat", "prob")
BGtest=TestBGpom
} # kon BG

# szereg to nazwa szeregu w ktorym badamy stalosc wariancji
# podzial to licznosc poszczegolnych podprobek np. podzial=c(20,25,30)
ConoverMulti=function(residuals, subsample)
{ # poczatek funkcji
szereg=residuals
podzial=subsample
# sprawdzenie poprawnosci podzialu
szereg=as.matrix(szereg)
licznosc=nrow(szereg)
PodzialSuma=sum(podzial)
if (PodzialSuma>licznosc)
 { # poczatek if
 podzial_pom=0
 for (i in 1:length(podzial))
  { # poczatek petli
  if (cumsum(podzial)[i]>licznosc)
  { # poczatek if
   if (sum(podzial_pom)<licznosc)
   { # poczatek if
   podzial_pom=podzial[1:(i-1)]
   podzial_pom[i]=licznosc-sum(podzial[1:(i-1)])
  } # koniec if
  } # koniec if
 } # koniec petli
 podzial=podzial_pom
 } # koniec if
if (PodzialSuma<licznosc)
 { # poczatek if
 podzial[length(podzial)]=podzial[length(podzial)]+licznosc-PodzialSuma
 } # koniec if
# koniec sprawdzania poprawnosci podzialu

# podzal na podszeregi
SzeregPodzielony=matrix(NA, max(podzial), length(podzial))
for (i in 1:ncol(SzeregPodzielony))
{ # poczatek petli
if (i==1)
 { ip=1
   ik=podzial[1] }
if (i!=1)
 { ip=cumsum(podzial)[i-1]+1
   ik=cumsum(podzial)[i] }
SzeregPodzielony_pom=as.matrix(szereg[ip:ik,])
if (nrow(SzeregPodzielony_pom)<nrow(SzeregPodzielony)) SzeregPodzielony_pom=rbind(SzeregPodzielony_pom, matrix(NA, (nrow(SzeregPodzielony)-nrow(SzeregPodzielony_pom)),1))
SzeregPodzielony[,i]=SzeregPodzielony_pom
} # koniec petli
# koniec podzialu na podszeregi

# liczenie odchylen od srednich
SrednieKolumn=matrix(NA,1,ncol(SzeregPodzielony))
SzeregPodzielonyOdchylenia=matrix(NA, max(podzial), length(podzial))
for (i in 1:ncol(SzeregPodzielony))
{ # poczatek petli
SrednieKolumn_pom=na.omit(SzeregPodzielony[,i])
SrednieKolumn[1,i]=sum(SrednieKolumn_pom)/length(SrednieKolumn_pom)    # srednia arytmetyczna danych w danej kolumnie
SzeregPodzielonyOdchylenia[,i]=SzeregPodzielony[,i]-SrednieKolumn[1,i] # odchylenia od sredniej
SzeregPodzielonyOdchylenia[,i]=abs(SzeregPodzielonyOdchylenia[,i])     # wartosc bezwzgledna odchylen
} # koniec petli

# tworzenie stakowanego szeregu wartosci bezwzglenych z odchylen
for (i in 1:ncol(SzeregPodzielony))
{ # poczatek petli
SzeregStackOdchylen_pom=as.matrix(na.omit(SzeregPodzielonyOdchylenia[,i]))
SzeregStackOdchylen_pom=cbind(SzeregStackOdchylen_pom, matrix(i, nrow(SzeregStackOdchylen_pom), ncol(SzeregStackOdchylen_pom)))
if (i==1) SzeregStackOdchylen=SzeregStackOdchylen_pom
if (i!=1) SzeregStackOdchylen=rbind(SzeregStackOdchylen, SzeregStackOdchylen_pom)
} # koniec petli

# sortowanie szeregu
kolejnosc=order(SzeregStackOdchylen[,1])
SzeregStackOdchylenSort=SzeregStackOdchylen[kolejnosc,]

# rangowanie
Rangi=rank(SzeregStackOdchylenSort[,1], ties.method="average")

# macierze: odchylenia, numer szeregu, rangi, rangi^2, rangi^4
Rangi=cbind(SzeregStackOdchylenSort[,1], SzeregStackOdchylenSort[,2], as.matrix(Rangi), as.matrix(Rangi^2), as.matrix(Rangi^4))

# liczenie si
Si=matrix(NA,1,ncol(SzeregPodzielony))
for (i in 1:ncol(SzeregPodzielony))
{ # poczatek petli
Si_pom<-Rangi[Rangi[,2]==i]
Si_pom=matrix(Si_pom, podzial[i], 5)
Si[1,i]=sum(Si_pom[,4])
} # koniec petli

# S_all
S_all=sum(Si)

# Sp
Sp=Si^2
Sp=Sp/t(as.matrix(podzial))
Sp=sum(Sp)

# Sr
Sr=sum(Rangi[,5])

# C
C=sum(Rangi[,4])^2/max(Rangi[,3])

# T
T=(max(Rangi[,3])-1)*(Sp-C)/(Sr-C)

# prob
lss=length(podzial)-1
prob=1-pchisq(T, lss)

ConoverMulti=t(as.matrix(c(T, prob)))
ConoverMulti=data.frame(ConoverMulti, row.names=c("Chi^2 stat"))
names(ConoverMulti)=c("stat", "prob")
ConoverMulti=data.frame(ConoverMulti, row.names=c("Chi^2 stat"))
} # koniec funkcji

JBtest=function(residuals)
{ # pocz f
szereg=residuals
m2=sum((szereg-mean(szereg))^2)/length(szereg)
pm2=sqrt(m2)
m3=sum((szereg-mean(szereg))^3)/length(szereg)
m4=sum((szereg-mean(szereg))^4)/length(szereg)
mi3=m3/pm2^3
mi4=m4/pm2^4
jb=length(szereg)*(1/6*mi3^2+1/24*(mi4-3)^2)
prob=1-pchisq(q=jb,df=2)
pomJBtest=data.frame(rbind(jb, prob), row.names=c("stat Chi^2","prob"))
names(pomJBtest)=""
JBtest=pomJBtest
} # kon f


GQtest=function(residuals, subsample, nep)
{ # pocz
szereg=residuals
podzial=subsample
lsp=nep
pp1=szereg[1:podzial[1]]
pp2=szereg[podzial[2]:length(szereg)]
w1=sum((pp1-mean(pp1))^2)/length(pp1)
w2=sum((pp2-mean(pp2))^2)/length(pp2)
if (w1>=w2)
{ gq=w1/w2
 lssl=length(pp1)
 lssm=length(pp2)}
if (w2>w1)
{ gq=w2/w1
 lssl=length(pp2)
 lssm=length(pp1)}
prob=1-pf(q=gq,df1=lssl-lsp, df2=lssm-lsp)
pomGQtest=data.frame(rbind(gq, prob), row.names=c("stat F","prob"))
names(pomGQtest)=""
GQtest=pomGQtest
} # kon

# paramTeta poczatkowa ocena parametrow dlugookresowych modelu (7)
# dateset macierz zawierajaca wszystkie zmienne potrzebne do oszacowania modelu
# vecSR wektor z numerami zmiennych wchodzacyhc w skalad relacji krotkookresowej, pierwszy MUSI byc numer (nazwa) zmiennej dy
# vecLR wektor z numerami zmiennych wchodzacyhc w skalad relacji dlugookresowej, pierwszy numer MUSI wskazywac ly
# quantities wektor licznosci T1,...,Tn poszczegolnych jednostek
# const=1 domyslnie, oznacza dolaczenie stalej do relacji SR
PMG=function(paramTeta, vecSR, vecLR, dataset, quantity, const)
{ # pocz bs
#paramTeta=c(1,-3,-0.5)
#vecSR=list(SR1=c(2,6,8), SR2=c(2,6,8),SR3=c(2,6,8),SR4=c(2,6,8),SR5=c(2,6,8),SR6=c(2,6,8),SR7=c(2,6,8),SR8=c(2,6,8),SR9=c(2,6,8))
#vecSR=list(SR1=c("dy10","dopeness","diip","ldcrisk","ddcpi"),SR2=c("dy10","dopeness","diip","ldcrisk","ddcpi"),SR3=c("dy10","dopeness","diip","ldcrisk","ddcpi"),SR4=c("dy10","dopeness","diip","ldcrisk","ddcpi"),SR5=c("dy10","dopeness","diip","ldcrisk","ddcpi"),SR6=c("dy10","dopeness","diip","ldcrisk","ddcpi"),SR7=c("dy10","dopeness","diip","ldcrisk","ddcpi"),SR8=c("dy10","dopeness","diip","ldcrisk","ddcpi"),SR9=c("dy10","dopeness","diip","ldcrisk","ddcpi"))
#vecLR=c(3,4,7,9)   # wektor wspolny dla wszystkich N
#vecLR=c("ly10","openess","iip","crisk","cpi")
#quantities=licznosci
#dataset=daneP
#const=TRUE
quantities=quantity
noN=length(quantities)

# rozpoznanie, czy vecSR i vacLR podane sa w postaci wektorow numerow zmiennych czy nazw zmiennych
# by program dziala poprawnie musza byc numery, ale chyba wygodniej zbiory danych tworzyc
# stosujac nazwy, dlatego musi byc mozliwosc konwersji z nazw na numery
if (!is.numeric(vecSR[[1]]))
{ # pocz if
pomvecSR=list()
for (i in 1:noN) { # pocz i
pomvecSR[[length(pomvecSR)+1]]=match(vecSR[[i]],names(dataset))
} # kon i
vecSR=pomvecSR
} # kon if
# teraz to samo z nazami w LR
if (!is.numeric(vecLR[1]))
{ # pocz if
pomvecLR=match(vecLR,names(dataset))
vecLR=pomvecLR
} # kon if

# zapamietanie pozycji zmiennej objasnianej
# usuniecie z wektora SR pozycji zmiennej objasnianej
nodY=NULL
vecSR2=vector("list", length=noN)
for (i in 1:noN)
{ # pocz i
nodY=c(nodY, unlist(vecSR[i])[1])
vecSR2[[i]]=as.numeric(unlist(vecSR[i])[-1])
} # kon i
# wektor paramTeta jako macierz kolumnowa
Teta=as.matrix(paramTeta)

# liczenie ksi teta relacji (7)
# liczenie dla kazdego i osobno, bo T1 nie musi byc rowne T2...Tn
# wektory SR tez nie musza byc takie same
# adresy poczatkow i koncow wierszy zmiennych jednostek panelu
do=cumsum(quantities)
od=(do-quantities)+1
KsiTeta=vector("list", length=noN)
ly=vector("list", length=noN)
X=vector("list", length=noN)
for (i in 1:noN)
{ # pocz i
# pobranie danych
pomly=as.matrix(dataset[od[i]:do[i],vecLR[1]])  # pierwsza zmienna na liscie ly
pomX=as.matrix(dataset[od[i]:do[i],vecLR[-1]])  # wszystkie LR wylaczywszy ly
# obliczenie KsiTeta (wzor 7)
pomKsiTeta=as.matrix(pomly-pomX%*%Teta)
# przekazanie wynikow do obiektu typu lista
ly[[i]]=pomly
X[[i]]=pomX
KsiTeta[[i]]=pomKsiTeta
} # kon i

# macierz (lista macierzy) W zmiennych relacji krotkookresowej
# macierz (lista macierzy) dy zmiennej objasnianej dy
# macierz (lista macierzy) wartosci teoretycznych
W=vector("list", length=noN)
dy=vector("list", length=noN)
yhat=vector("list", length=noN)
for (i in 1:noN)
{ # pocz i
pomW=as.matrix(dataset[od[i]:do[i],vecSR2[[i]]])
if (const==1) pomW=cbind(pomW,1)
W[[i]]=pomW
dy[[i]]=as.matrix(dataset[od[i]:do[i],nodY[i]])
} # kon i

# macierz H=I-W(W'W)^-1W'
H=vector("list", length=noN)

# parametry ecm wedlug relacji (10)
Fi=vector("list", length=noN)

# wariancje sigma^2 wedlug relacji (11)
sigma=vector("list", length=noN)

# nawias1 i nawias2 relacji (9)
Nawias1=vector("list", length=noN)
Nawias2=vector("list", length=noN)

# sumy do LogL
suma1=0
suma2=0

# parametry i bledy relacji krotkookresowej
paramSR=vector("list", length=noN)

# reszty relacji (6)
reszty=vector("list", length=noN)

# obiekty do estymatora macierz wariancji i kowariancji czyli elementy diagonalne
# i piewszy wiersz zkladajacy sie z varTeta i varFirstRow (por strona 626)
varTeta=matrix(0, (length(vecLR)-1), (length(vecLR)-1))
varFirstRowFi=NULL
varFirstRowSR=NULL
varFi=matrix(0, noN, noN)
varSR=NULL
varMixedFiSR=NULL

# adresy koncow varFirstRowSR w macierzy  macierzy wariancji i kowariancji
# one zaleza od liczby kolumn (zmiennych) w macierzy W ktore maga sie roznic
# miedzy jednostkami panelu
doKolSR=NULL

#i=1
for (i in 1:noN)
{ # pocz i
pomW=W[[i]]                                                                         # tworzenie macierzy W
pomH=diag(1,nrow(pomW), nrow(pomW))-pomW%*%solve(t(pomW)%*%pomW)%*%t(pomW)          # tworzeniem macierzy H
H[[i]]=pomH                                                                         # wybrane Hi

pomKsiTeta=KsiTeta[[i]]
pomdy=dy[[i]]

pomFi=solve(t(pomKsiTeta)%*%pomH%*%pomKsiTeta)%*%t(pomKsiTeta)%*%pomH%*%pomdy       # parametry fi liczone z relacji (10)
pomFi=as.numeric(pomFi)
Fi[[i]]=pomFi

pomSigma=t((pomdy-pomKsiTeta*pomFi))%*%pomH%*%(pomdy-pomKsiTeta*pomFi)/nrow(pomW)   # sigma^2 liczona z relacji (11)
sigma[[i]]=pomSigma

pomX=X[[i]]
pomNawias1=(t(pomX)%*%pomH%*%pomX)*as.numeric((pomFi^2/pomSigma))                  # skladowa i nawiasu pierwszegro relacji (9)
Nawias1[[i]]=pomNawias1

pomly=ly[[i]]
pomNawias2=as.numeric((pomFi/pomSigma))*t(pomX)%*%pomH%*%(pomdy-pomFi*pomly)        # skladowa i nawiasu drugiego relacji (9)
Nawias2[[i]]=pomNawias2

# liczenie LogL
suma1=suma1+log(2*pi*pomSigma)
suma2=suma2+(1/pomSigma)*t((pomdy-pomFi*pomKsiTeta))%*%pomH%*%(pomdy-pomFi*pomKsiTeta)

# estymacja parametrow czesci SR
pomModel=lm(pomdy~pomKsiTeta+pomW-1)       # pomW ma juz kolumne jedynek, dlatego -1
pomSR=as.matrix(coef(pomModel)[-1])        # usuniecie ponownego oszacowania Fi
paramSR[[i]]=pomSR

# obliczanie reszt
# najpierw yhat
pomyhat=pomFi*pomKsiTeta+pomW%*%pomSR
pomReszta=as.matrix(pomdy-pomyhat)
reszty[[i]]=pomReszta
yhat[[i]]=pomyhat

# estymacja diagonalnej macierzy wariancji i kowariancji
varTeta=varTeta+as.numeric((pomFi^2/pomSigma))*(t(pomX)%*%pomX)
varFirstRowFi=cbind(varFirstRowFi,(t(pomX)%*%pomKsiTeta)*as.numeric(pomFi/pomSigma))
varFirstRowSR=cbind(varFirstRowSR,(t(pomX)%*%pomW)*as.numeric(pomFi/pomSigma))
doKolSR=c(doKolSR,ncol(varFirstRowSR))
varFi[i,i]=(t(pomKsiTeta)%*%pomKsiTeta)*as.numeric(1/pomSigma)
if (i==1)
{ # pocz if
varSR=(t(pomW)%*%pomW)*as.numeric(1/pomSigma)
varMixedFiSR=(t(pomKsiTeta)%*%pomW)*as.numeric(1/pomSigma)
} else
{ # pocz else
pomVarSR=(t(pomW)%*%pomW)*as.numeric(1/pomSigma)
pomVarMixedFiSR=(t(pomKsiTeta)%*%pomW)*as.numeric(1/pomSigma)
# dolaczanie macierzy zerowych do varSR
pom01VarSR=matrix(0,nrow(varSR),ncol(pomVarSR))  # dodana z prawej
pom02VarSR=matrix(0,nrow(pomVarSR),ncol(varSR))  # dodana z dolu
varSR=cbind(varSR,pom01VarSR)
varSR=rbind(varSR,cbind(pom02VarSR,pomVarSR))
# dolaczenie macierzy zerowych do varMixedFiSR
pom01VarMixedFiSR=matrix(0,nrow(varMixedFiSR),ncol(pomVarMixedFiSR))  # dodana z prawej
pom02VarMixedFiSR=matrix(0,1,ncol(varMixedFiSR))  # dodana z dolu
varMixedFiSR=cbind(varMixedFiSR,pom01VarMixedFiSR)
varMixedFiSR=rbind(varMixedFiSR,cbind(pom02VarMixedFiSR,pomVarMixedFiSR))
} # kon else
} # kon i

# nowe oceny parametrow Teta
pomNawias1=matrix(0,nrow(Teta),nrow(Teta))
pomNawias2=matrix(0,nrow(Teta),1)
for (i in 1:noN)
{ # pocz i
pomNawias1=pomNawias1+as.matrix(Nawias1[[i]])
pomNawias2=pomNawias2+as.matrix(Nawias2[[i]])
} # kon i
TetaHat=(-1)*solve(pomNawias1)%*%pomNawias2   # nowa ocena parametru teta z relacji (9)

# liczenie LogL
LogL=(-nrow(pomW)/2)*suma1-0.5*suma2

# liczenie macierzy wariancji i kowariancji
wiersz1=cbind(varFirstRowFi, varFirstRowSR)  # bez elementu 1,1 czyli tety
wiersz2=cbind(varFi, varMixedFiSR)
wiersz3=cbind(t(varMixedFiSR), varSR)
kolumna1=t(cbind(varTeta, varFirstRowFi, varFirstRowSR))
varcovarAll=cbind(kolumna1, rbind(wiersz1, wiersz2, wiersz3))
varcovarAllodw=solve(varcovarAll)
BledySzacunku=sqrt(diag(varcovarAllodw))

# liczenie statystyk t-studenta i prob w czesci dlugookresowej
paramLR=TetaHat
bledyLR=as.matrix(BledySzacunku[1:nrow(TetaHat)])
tstatLR=as.matrix(paramLR/bledyLR)
probLR=as.matrix((1-pnorm(abs(tstatLR), mean=0, sd=1, lower.tail=TRUE))*2)
Low95LR=as.matrix(paramLR+qnorm(p=0.025, mean=0, sd=1)*bledyLR)
High95LR=as.matrix(paramLR-qnorm(p=0.025, mean=0, sd=1)*bledyLR)
paramLR=data.frame(cbind(paramLR,bledyLR, tstatLR, probLR, Low95LR, High95LR))
names(paramLR)=c("Coef","StdErr","z","P>|z|","Low95%","High95%")

# dodawanie nazw do obiektow fi i sigma
Fi2=as.matrix(unlist(Fi))
sigma2=as.matrix(unlist(sigma))
row.names(Fi2)=paste(rep("ec",nrow(Fi2)), seq(from=1, to=nrow(Fi2), by=1), sep="")
row.names(sigma2)=paste(rep("sigma2_",nrow(sigma2)), seq(from=1, to=nrow(sigma2), by=1), sep="")
varianceLR=data.frame(sigma2)
names(varianceLR)="varianceLR"
# dodawanie nazw do zmiennych relacji krotkookresowej
VariableNames=names(dataset)
if (is.null(VariableNames)) VariableNames=paste(rep("VariableName",ncol(dataset)),1:ncol(dataset),sep="")
SRNames=vector("list", length=noN)
for (i in 1:noN)
{ # pocz i
pomSRNames=VariableNames[vecSR2[[i]]]
if (const) pomSRNames=c(pomSRNames,"const")
SRNames[[i]]=pomSRNames
row.names(paramSR[[i]])=pomSRNames
} # kon i

# liczenie statystyk t-studenta i prob w czesci krotkookresowej
# wyznaczanie adresow czesci fi i czesci SR w obiekcie BledySzacunku
# gdzie zaczynaja sie bledy fi i gdzie sa ich numery na przekatenej macierzy
# wariancji i kowariancji i tym sammym obiektu BledySzacunku
oddoFi=(nrow(TetaHat)+1):(nrow(TetaHat)+noN)
# adresy dla SR
oddoSR=list()
oddoSR[[length(oddoSR)+1]]<-seq(from=max(oddoFi)+1, len=nrow(paramSR[[1]]), by=1)
for (i in 2:noN)
{ # pocz i
oddoSR[[length(oddoSR)+1]]<-seq(from=max(oddoSR[[i-1]])+1, len=nrow(paramSR[[i]]), by=1)
} # kon i
paramSR2=list()
for (i in 1:noN)
{ # pocz i
pom_paramSR=rbind(as.matrix(Fi2[i,1]),as.matrix(paramSR[[i]]))
pom_bledySR=rbind(as.matrix(BledySzacunku[oddoFi[i]]),as.matrix(BledySzacunku[oddoSR[[i]]]))
pom_tstatSR=as.matrix(pom_paramSR/pom_bledySR)
probSR=as.matrix((1-pnorm(abs(pom_tstatSR), mean=0, sd=1, lower.tail=TRUE))*2)
Low95SR=as.matrix(pom_paramSR+qnorm(p=0.025, mean=0, sd=1)*pom_bledySR)
High95SR=as.matrix(pom_paramSR-qnorm(p=0.025, mean=0, sd=1)*pom_bledySR)
pom_paramSR=data.frame(cbind(pom_paramSR,pom_bledySR, pom_tstatSR, probSR, Low95SR, High95SR))
names(pom_paramSR)=c(paste("Coef(i=",i,")", sep=""),"StdErr","z","P>|z|","Low95%","High95%")
paramSR2[[length(paramSR2)+1]]<-pom_paramSR
} # kon i
paramSR=paramSR2

# liczenie reszt i testy reszt
TestyDiag=matrix(NA, noN, 12)
for (i in 1:noN)
{ # pocz i
reszta01=reszty[[i]]
podzialGQ=round(nrow(reszta01)/2,0)
podzialGQ=c(podzialGQ, (podzialGQ+1))                             # podzial na podproby w tescie GQ, numer ostatni pierwszej podproby i numer pierwszy drugiej podproby
podzialCON=round(nrow(reszta01)/3,0)
podzialCON=c(podzialCON,podzialCON,(nrow(reszta01)-2*podzialCON)) # podzial dla conovera, liczebnosci podprob
lsp=nrow(paramSR[[i]])                                            # liczba oszacowanych parametrow wraz z teta
X_testBG=cbind(KsiTeta[[i]],W[[i]])                               # macierz X do testu Breuscha-Godfreya, z wylaczeniem kolumny jedynek
if (const) X_testBG=as.matrix(X_testBG[,-ncol(X_testBG)])
acor.ord=4                                                        # max testowana autokorelacja w BG
pomX_toBG=as.matrix(X_testBG)
wJB=JBtest(residuals=reszta01)
TestyDiag[i,3]=wJB[1,1]
TestyDiag[i,4]=wJB[2,1]
wGQ=GQtest(residuals=reszta01, subsample=podzialGQ, nep=lsp)
TestyDiag[i,5]=wGQ[1,1]
TestyDiag[i,6]=wGQ[2,1]
wCON=ConoverMulti(residuals=reszta01, subsample=podzialCON)
TestyDiag[i,7]=wCON[1,1]
TestyDiag[i,8]=wCON[1,2]
pomdy=as.matrix(dy[[i]])
pomyhat=as.matrix(yhat[[i]])
R2=(cor(pomdy, pomyhat, method="pearson"))^2
R2bar=1-((1-R2)*(nrow(reszta01)-1)/(nrow(reszta01)-lsp))
Se=sqrt(sum((reszta01)^2)/(length(reszta01)-lsp))
V=Se/mean(pomdy)
TestyDiag[i,1]=R2
TestyDiag[i,2]=R2bar
wBG=BGtest(residuals=reszta01, explvariab=pomX_toBG, acor.ord=acor.ord)
TestyDiag[i,9]=wBG[1,1]         # test autokorelacji BG wersja F statysytyka
TestyDiag[i,10]=wBG[1,2]        # test autokorelacji BG wersja F prob
TestyDiag[i,11]=wBG[2,1]        # test autokorelacji BG wersja Chi^2 statysytyka
TestyDiag[i,12]=wBG[2,2]        # test autokorelacji BG wersja Chi^2 prob
} # kon i
NazwyWierszy=paste(rep("i=",noN), seq(from=1, to=noN, by=1), sep="")
TestyDiag=data.frame(TestyDiag, row.names=NazwyWierszy)
names(TestyDiag)=c("R^2","Rbar^2","statJB", "probJB", "statGQ", "probGQ", "statCON", "probCON",
"statBG(F)", "probBG(F)", "statBG(Chi^2)", "probBG(Chi^2)")
TestyDiag
PMG=list(TetaHat=TetaHat, LogL=LogL, SR=paramSR, LR=paramLR, varianceLR=varianceLR, DiagTests=TestyDiag, residuals=reszty)
} # kon bs


optimPMG=function(dLL, maxIter, TetaStart, vecSR, vecLR, dataset, quantity, const)
{ # pocz pmg
i=0
dLLHat=dLL+1
while (i<=maxIter & dLLHat>=dLL)
{ # pocz while
if (i==0) poBS1=PMG(paramTeta=TetaStart, vecSR=vecSR, vecLR=vecLR, dataset=dataset, quantity=quantity, const=const)
if (i>0) poBS1=PMG(paramTeta=poBS2$TetaHat, vecSR=vecSR, vecLR=vecLR, dataset=dataset,quantity=quantity, const=const)
poLL1=poBS1$LogL
poBS2=PMG(paramTeta=poBS1$TetaHat, vecSR=vecSR, vecLR=vecLR, dataset=dataset,quantity=quantity, const=const)
poLL2=poBS2$LogL
dLLHat=poLL2-poLL1
i=i+1
} # kon while
optimPMG=list(LogL=poLL2, dLogL=dLLHat, i=i, LR=poBS2$LR, SR=poBS2$SR, DiagTests=poBS2$DiagTests, residuals=poBS2$residuals)
} # kon pmg

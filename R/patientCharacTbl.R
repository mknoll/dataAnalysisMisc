#' @title Creates a patient characteristics table
#' 
#' @import plyr
crtPatTbl <- function(data) {
  out <- list()
  i<-1
  ## Alter-Kategorien
  if ("AGE" %in% colnames(data)) {
    ageCat <- c(0, 50,60,70, 200)
    ages <- list()
    u <- 1
    for (j in 2:(length(ageCat)-1)) {
      vec <- data.frame(name=paste("<=",ageCat[j],sep=""), 
                        n=length(data$AGE[which(data$AGE >= ageCat[j-1] & data$AGE <= ageCat[j])]))
      ages[[u]] <- vec
      u<-u+1
    }
    ages[[u]] <- data.frame(name=paste(">",ageCat[j],sep=""), 
                      n=length(data$AGE[which(data$AGE > ageCat[j])]))
    
    out[[i]] <- list(name="AGE", data=ages)
    i<-i+1
  }
  ##Geschlecht
  if ("GENDER" %in% colnames(data)) {
    gen <- data.frame(name=c("M", "F"), n=c(0,0))
    gen[which(gen$name == "M"),"n"] <- length(data[which(data$GENDER == "M"),1])
    gen[which(gen$name == "F"),"n"] <- length(data[which(data$GENDER == "F"),1])
    if (length(data[,1])-gen$n[1]-gen$n[2] != 0) {
      gen <- rbind.fill(name=gen, n=data.frame("NA", abs(length(data[,1])-gen$n[1]-gen$n[2])))
    }
    out[[i]] <- list(name="GENDER", data=gen)
    i <- i+1
  }
  ##KPS
  if ("KPS" %in% colnames(data)) {
    kps <- data.frame(name=c("<=80", ">80"), n=c(0,0)) 
    kps[which(kps$name == "<=80"),"n"] <- length(data[which(data$KPS <= 80),1])
    kps[which(kps$name == ">80"),"n"] <- length(data[which(data$KPS  > 80),1])
    if (length(data[,1])-kps$n[1]-kps$n[2] != 0) {
      kps <- rbind.fill(kps, data.frame(name="NA", n=abs(length(data[,1])-kps$n[1]-kps$n[2])))
    }
    out[[i]] <- list(name="KPS", data=kps)
    i <- i+1
  }
  ##IDH1 mut
  if ("IDH1" %in% colnames(data)) {
    idh <- data.frame(name=c("wt", "mut"), n=c(0,0))
    idh[which(idh$name == "wt"),"n"] <- length(data[which(data$IDH1 == "wt"),1])
    idh[which(idh$name == "mut"),"n"] <- length(data[which(data$IDH1 == "mut" ),1])
    if (length(data[,1])-idh$n[1]-idh$n[2] != 0) {
      idh <- rbind.fill(idh, data.frame(name="NA", n=abs(length(data[,1])-idh$n[1]-idh$n[2])))
    }
    out[[i]] <- list(name="IDH1", data=idh)
    i <- i+1
  }
  ##MGMT Hypermethylation
  if ("MGMT" %in% colnames(data)) {
    mgmt <- data.frame(name=c("hypermeth", "non-hypermeth"), n=c(0,0))
    mgmt[which(mgmt$name == "hypermeth"),"n"] <- length(data[which(data$MGMT == 1),1])
    mgmt[which(mgmt$name == "non-hypermeth"),"n"] <- length(data[which(data$MGMT == 0),1])
    if (length(data[,1])-length(data[which(data$MGMT == 1),1])-length(data[which(data$MGMT == 0),1]) != 0) {
      mgmt <- rbind.fill(mgmt, data.frame(name="NA", n=abs(length(data[,1])-mgmt$n[1]-mgmt$n[2])))
    }
    out[[i]] <- list(name="MGMT", data=mgmt)
    i <- i+1
  }
  ##Gesamtdosis
  if ("RT_GESAMT" %in% colnames(data)) {
    rt <- data.frame(name=c("Dose"), 
                     lage=mean(data$RT_GESAMT, na.rm=T), 
                     streuung=sd(data$RT_GESAMT, na.rm=T), 
                     na=length(data$RT_GESAMT[is.na(data$RT_GESAMT)]))
    out[[i]] <- list(name="RT total", data=rt)
    i <- i+1
  }
  ##Anzahl Fraktionen
  if ("FX" %in% colnames(data)) {
    fx <- data.frame(name=c("Fx"), 
                     lage=modus(data$FX), 
                     streuung=c(min(data$FX, na.rm=T), max(data$FX, na.rm=T)), 
                     na=length(data$FX[is.na(data$FX)]))
    out[[i]] <- list(name="FX", data=fx)
    i <- i+1
  }
  ##GTV Vol
  if ("GTV" %in% colnames(data)) {
    gtv <- data.frame(name=c("GTV"), 
                     lage=mean(data$GTV, na.rm=T), 
                     streuung=sd(data$GTV, na.rm=T), 
                     na=length(data$GTV[is.na(data$GTV)]))
    out[[i]] <- list(name="GTV Volume", data=gtv)
    i <- i+1
  }
  ##TUR
  if ("TUR" %in% colnames(data)) {
    tur <- data.frame(name=c("TUR"), 
                      lage=mean(data$TUR, na.rm=T), 
                      streuung=sd(data$TUR, na.rm=T), 
                      na=length(data$TUR[is.na(data$TUR)]))
    out[[i]] <- list(name="TUR", data=tur)
    i <- i+1
  }
  ## Total
  out[[i]] <- list(name="Total", data=data.frame(name="Total", n=length(data[,1])))
  i<-i+1
  
  return(out)
}



modus <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


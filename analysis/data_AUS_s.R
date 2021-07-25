#setwd("/am/miro/home/iliu/TEX/paper/DIF/PISA2012")
#
# Upload Sam's csv data format
#
pisa.AUS<-read.csv(file("australian_student_science_questions.csv"))
pisa.AUS[1:10,]
summary(pisa.AUS)
nrow(pisa.AUS)
pisa.AUS.s<-pisa.AUS[which(pisa.AUS$CognitiveItemCode=='PS131Q02D' |  # Only keep these 53 questions
                           pisa.AUS$CognitiveItemCode=='PS131Q04D' |
                           pisa.AUS$CognitiveItemCode=='PS256Q01' |
                           pisa.AUS$CognitiveItemCode=='PS269Q01' |
                           pisa.AUS$CognitiveItemCode=='PS269Q03D' |
                           pisa.AUS$CognitiveItemCode=='PS269Q04T' |
                           pisa.AUS$CognitiveItemCode=='PS326Q01' |
                           pisa.AUS$CognitiveItemCode=='PS326Q02' |
                           pisa.AUS$CognitiveItemCode=='PS326Q03' |
                           pisa.AUS$CognitiveItemCode=='PS326Q04T' |
                           pisa.AUS$CognitiveItemCode=='PS408Q01' |
                           pisa.AUS$CognitiveItemCode=='PS408Q03' |
                           pisa.AUS$CognitiveItemCode=='PS408Q04T' |
                           pisa.AUS$CognitiveItemCode=='PS408Q05' |
                           pisa.AUS$CognitiveItemCode=='PS413Q04T' |
                           pisa.AUS$CognitiveItemCode=='PS413Q05' |
                           pisa.AUS$CognitiveItemCode=='PS413Q06' |
                           pisa.AUS$CognitiveItemCode=='PS415Q02' |
                           pisa.AUS$CognitiveItemCode=='PS415Q07T' |
                           pisa.AUS$CognitiveItemCode=='PS415Q08T' |
                           pisa.AUS$CognitiveItemCode=='PS425Q02' |
                           pisa.AUS$CognitiveItemCode=='PS425Q03' |
                           pisa.AUS$CognitiveItemCode=='PS425Q04' |
                           pisa.AUS$CognitiveItemCode=='PS425Q05' |
                           pisa.AUS$CognitiveItemCode=='PS428Q01' |
                           pisa.AUS$CognitiveItemCode=='PS428Q03' |
                           pisa.AUS$CognitiveItemCode=='PS428Q05' |
                           pisa.AUS$CognitiveItemCode=='PS438Q01T' |
                           pisa.AUS$CognitiveItemCode=='PS438Q02' |
                           pisa.AUS$CognitiveItemCode=='PS438Q03D' |
                           pisa.AUS$CognitiveItemCode=='PS465Q01' |
                           pisa.AUS$CognitiveItemCode=='PS465Q02' |
                           pisa.AUS$CognitiveItemCode=='PS465Q04' |
                           pisa.AUS$CognitiveItemCode=='PS466Q01T' |
                           pisa.AUS$CognitiveItemCode=='PS466Q05' |
                           pisa.AUS$CognitiveItemCode=='PS466Q07T' |
                           pisa.AUS$CognitiveItemCode=='PS478Q01' |
                           pisa.AUS$CognitiveItemCode=='PS478Q02T' |
                           pisa.AUS$CognitiveItemCode=='PS478Q03T' |
                           pisa.AUS$CognitiveItemCode=='PS498Q02T' |
                           pisa.AUS$CognitiveItemCode=='PS498Q03' |
                           pisa.AUS$CognitiveItemCode=='PS498Q04' |
                           pisa.AUS$CognitiveItemCode=='PS514Q02' |
                           pisa.AUS$CognitiveItemCode=='PS514Q03' |
                           pisa.AUS$CognitiveItemCode=='PS514Q04' |
                           pisa.AUS$CognitiveItemCode=='PS519Q01' |
                           pisa.AUS$CognitiveItemCode=='PS519Q02T' |
                           pisa.AUS$CognitiveItemCode=='PS519Q03' |
                           pisa.AUS$CognitiveItemCode=='PS521Q02' |
                           pisa.AUS$CognitiveItemCode=='PS521Q06' |
                           pisa.AUS$CognitiveItemCode=='PS527Q01T' |
                           pisa.AUS$CognitiveItemCode=='PS527Q03T' |
                           pisa.AUS$CognitiveItemCode=='PS527Q04T') ,]
pisa.AUS.s[1:10,]
nrow(pisa.AUS.s)/53  #number of students in the data
#
library(reshape)
library(data.table)
library(tidyr)
setDT(pisa.AUS.s)
pisa.AUS.s.w<-dcast(pisa.AUS.s, StIDStd + Gender ~ CognitiveItemCode, value.var="CognitiveItemScoreBoolean") #wide format and keep gender and ID
pisa.AUS.s.w[1:10,]
nrow(pisa.AUS.s.w)
summary(pisa.AUS.s.w)
#
row.has.na<-apply(pisa.AUS.s.w,1,function(x){any(is.na(x))})
sum(row.has.na)  # every students have at lesat one NA
final.data<-pisa.AUS.s.w[!row.has.na,]  # no one left
#
pisa.AUS.s.w1<-dcast(pisa.AUS.s, StIDStd + Gender ~ CognitiveItemCode, value.var="CognitiveItemScore") #wide format and keep gender and ID
pisa.AUS.s.w1[1:10,]
nrow(pisa.AUS.s.w1)
row.has.seven<-apply(pisa.AUS.s.w1,1,function(x){any(x==7)})
sum(row.has.seven) # number of students that have at lesat one 7 response.
final.data.1<-pisa.AUS.s.w1[!row.has.seven,]
final.data.1 # every students have at lesat one 7 response.
summary(pisa.AUS.s.w1)

table(pisa.AUS.s$SCHOOLID)/53

#School 1
pisa.AUS.s1<-pisa.AUS.s[which(pisa.AUS.s$SCHOOLID=='1') ,]
pisa.AUS.s1[1:10,]
nrow(pisa.AUS.s1)/53
pisa.AUS.s1.w1<-dcast(pisa.AUS.s1, StIDStd + Gender ~ CognitiveItemCode, value.var="CognitiveItemScore") #wide format and keep gender and ID
pisa.AUS.s1.w1[1:10,]
nrow(pisa.AUS.s1.w1)


#BOOKID 1  - The students who took BOOKID=1
pisa.AUS.b1<-pisa.AUS.s[which(pisa.AUS.s$BOOKID=='1') ,]
pisa.AUS.b1[1:10,]
nrow(pisa.AUS.b1)/53
pisa.AUS.b1.w1<-dcast(pisa.AUS.b1, StIDStd + Gender ~ CognitiveItemCode, value.var="CognitiveItemScore") #wide format and keep gender and ID
pisa.AUS.b1.w1[1:10,]
nrow(pisa.AUS.b1.w1)
summary(pisa.AUS.b1.w1)
#    StIDStd         Gender            PS131Q02D        PS131Q04D         PS256Q01        PS269Q01   PS269Q03D   PS269Q04T
# Min.   :   11   Length:1145        Min.   :0.0000   Min.   :0.0000   Min.   :0.000   Min.   :7   Min.   :7   Min.   :7
# 1st Qu.: 3582   Class :character   1st Qu.:0.0000   1st Qu.:0.0000   1st Qu.:1.000   1st Qu.:7   1st Qu.:7   1st Qu.:7
# Median : 7123   Mode  :character   Median :1.0000   Median :0.0000   Median :1.000   Median :7   Median :7   Median :7
# Mean   : 7174                      Mean   :0.6131   Mean   :0.3956   Mean   :1.122   Mean   :7   Mean   :7   Mean   :7
# 3rd Qu.:10751                      3rd Qu.:1.0000   3rd Qu.:1.0000   3rd Qu.:1.000   3rd Qu.:7   3rd Qu.:7   3rd Qu.:7
# Max.   :14480                      Max.   :8.0000   Max.   :8.0000   Max.   :8.000   Max.   :7   Max.   :7   Max.   :7
#    PS326Q01         PS326Q02         PS326Q03        PS326Q04T         PS408Q01    PS408Q03   PS408Q04T    PS408Q05
# Min.   :0.0000   Min.   :0.0000   Min.   :0.0000   Min.   :0.0000   Min.   :7   Min.   :7   Min.   :7   Min.   :7
# 1st Qu.:0.0000   1st Qu.:0.0000   1st Qu.:0.0000   1st Qu.:0.0000   1st Qu.:7   1st Qu.:7   1st Qu.:7   1st Qu.:7
# Median :1.0000   Median :1.0000   Median :1.0000   Median :0.0000   Median :7   Median :7   Median :7   Median :7
# Mean   :0.8288   Mean   :0.8681   Mean   :0.7764   Mean   :0.5048   Mean   :7   Mean   :7   Mean   :7   Mean   :7
# 3rd Qu.:1.0000   3rd Qu.:1.0000   3rd Qu.:1.0000   3rd Qu.:1.0000   3rd Qu.:7   3rd Qu.:7   3rd Qu.:7   3rd Qu.:7
# Max.   :8.0000   Max.   :8.0000   Max.   :8.0000   Max.   :8.0000   Max.   :7   Max.   :7   Max.   :7   Max.   :7
#   PS413Q04T         PS413Q05         PS413Q06         PS415Q02        PS415Q07T       PS415Q08T         PS425Q02
# Min.   :0.0000   Min.   :0.0000   Min.   :0.0000   Min.   :0.0000   Min.   :0.000   Min.   :0.0000   Min.   :0.000
# 1st Qu.:0.0000   1st Qu.:0.0000   1st Qu.:0.0000   1st Qu.:1.0000   1st Qu.:1.000   1st Qu.:0.0000   1st Qu.:0.000
# Median :1.0000   Median :1.0000   Median :0.0000   Median :1.0000   Median :1.000   Median :1.0000   Median :1.000
# Mean   :0.8061   Mean   :0.9127   Mean   :0.6035   Mean   :0.9293   Mean   :0.959   Mean   :0.7205   Mean   :1.167
# 3rd Qu.:1.0000   3rd Qu.:1.0000   3rd Qu.:1.0000   3rd Qu.:1.0000   3rd Qu.:1.000   3rd Qu.:1.0000   3rd Qu.:1.000
# Max.   :8.0000   Max.   :8.0000   Max.   :8.0000   Max.   :8.0000   Max.   :8.000   Max.   :8.0000   Max.   :8.000
#    PS425Q03         PS425Q04        PS425Q05        PS428Q01         PS428Q03         PS428Q05        PS438Q01T
# Min.   :0.0000   Min.   :0.000   Min.   :0.000   Min.   :0.0000   Min.   :0.0000   Min.   :0.0000   Min.   :0.0000
# 1st Qu.:0.0000   1st Qu.:0.000   1st Qu.:0.000   1st Qu.:0.0000   1st Qu.:1.0000   1st Qu.:0.0000   1st Qu.:1.0000
# Median :0.0000   Median :0.000   Median :1.000   Median :1.0000   Median :1.0000   Median :0.0000   Median :1.0000
# Mean   :0.9293   Mean   :1.014   Mean   :1.271   Mean   :0.6376   Mean   :0.8908   Mean   :0.5651   Mean   :0.9039
# 3rd Qu.:1.0000   3rd Qu.:1.000   3rd Qu.:1.000   3rd Qu.:1.0000   3rd Qu.:1.0000   3rd Qu.:1.0000   3rd Qu.:1.0000
# Max.   :8.0000   Max.   :8.000   Max.   :8.000   Max.   :8.0000   Max.   :8.0000   Max.   :8.0000   Max.   :8.0000
#    PS438Q02        PS438Q03D         PS465Q01        PS465Q02         PS465Q04        PS466Q01T    PS466Q05   PS466Q07T
# Min.   :0.0000   Min.   :0.0000   Min.   :0.000   Min.   :0.0000   Min.   :0.0000   Min.   :7   Min.   :7   Min.   :7
# 1st Qu.:1.0000   1st Qu.:0.0000   1st Qu.:0.000   1st Qu.:0.0000   1st Qu.:0.0000   1st Qu.:7   1st Qu.:7   1st Qu.:7
# Median :1.0000   Median :1.0000   Median :1.000   Median :1.0000   Median :0.0000   Median :7   Median :7   Median :7
# Mean   :0.8882   Mean   :0.6192   Mean   :1.071   Mean   :0.7598   Mean   :0.5013   Mean   :7   Mean   :7   Mean   :7
# 3rd Qu.:1.0000   3rd Qu.:1.0000   3rd Qu.:2.000   3rd Qu.:1.0000   3rd Qu.:1.0000   3rd Qu.:7   3rd Qu.:7   3rd Qu.:7
# Max.   :8.0000   Max.   :8.0000   Max.   :8.000   Max.   :8.0000   Max.   :8.0000   Max.   :7   Max.   :7   Max.   :7
#    PS478Q01        PS478Q02T        PS478Q03T        PS498Q02T      PS498Q03         PS498Q04        PS514Q02
# Min.   :0.0000   Min.   :0.0000   Min.   :0.0000   Min.   :0.0   Min.   :0.0000   Min.   :0.000   Min.   :0.0000
# 1st Qu.:0.0000   1st Qu.:0.0000   1st Qu.:0.0000   1st Qu.:0.0   1st Qu.:0.0000   1st Qu.:0.000   1st Qu.:1.0000
# Median :1.0000   Median :1.0000   Median :1.0000   Median :0.0   Median :0.0000   Median :2.000   Median :1.0000
# Mean   :0.8341   Mean   :0.7537   Mean   :0.9065   Mean   :0.8   Mean   :0.7991   Mean   :1.672   Mean   :0.9607
# 3rd Qu.:1.0000   3rd Qu.:1.0000   3rd Qu.:1.0000   3rd Qu.:1.0   3rd Qu.:1.0000   3rd Qu.:2.000   3rd Qu.:1.0000
# Max.   :8.0000   Max.   :8.0000   Max.   :8.0000   Max.   :8.0   Max.   :8.0000   Max.   :8.000   Max.   :8.0000
#    PS514Q03         PS514Q04         PS519Q01   PS519Q02T    PS519Q03    PS521Q02    PS521Q06   PS527Q01T   PS527Q03T
# Min.   :0.0000   Min.   :0.0000   Min.   :7   Min.   :7   Min.   :7   Min.   :7   Min.   :7   Min.   :7   Min.   :7
# 1st Qu.:0.0000   1st Qu.:0.0000   1st Qu.:7   1st Qu.:7   1st Qu.:7   1st Qu.:7   1st Qu.:7   1st Qu.:7   1st Qu.:7
# Median :1.0000   Median :1.0000   Median :7   Median :7   Median :7   Median :7   Median :7   Median :7   Median :7
# Mean   :0.6341   Mean   :0.6358   Mean   :7   Mean   :7   Mean   :7   Mean   :7   Mean   :7   Mean   :7   Mean   :7
# 3rd Qu.:1.0000   3rd Qu.:1.0000   3rd Qu.:7   3rd Qu.:7   3rd Qu.:7   3rd Qu.:7   3rd Qu.:7   3rd Qu.:7   3rd Qu.:7
# Max.   :8.0000   Max.   :8.0000   Max.   :7   Max.   :7   Max.   :7   Max.   :7   Max.   :7   Max.   :7   Max.   :7
#   PS527Q04T
# Min.   :7
# 1st Qu.:7
# Median :7
# Mean   :7
# 3rd Qu.:7
# Max.   :7


pisa.AUS.book1<-pisa.AUS.b1.w1[,c(1:5,9:12,17:35,39:47)]
pisa.AUS.book1[1:10,]
summary(pisa.AUS.book1)
ncol(pisa.AUS.book1)
#37 (35 science items)
nrow(pisa.AUS.book1)
#1145 students

save(pisa.AUS.book1, file="pisa_AUS_bookID_1.RData")






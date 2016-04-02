
## households
ss13husa = read.csv("~/Documents/Altman_HD/2013-american-community-survey/pums/ss13husa.csv")

for(a in 1:231){
    print(class(ss13husa[,a]))
}

## v large
ss13pusa = read.csv("~/Documents/Altman_HD/2013-american-community-survey/pums/ss13pusa.csv")



Runyear = 2016

if Runyear == 2016:
    import Samplelist as Slist
elif Runyear == 2017:
    #import Samplelist2017 as Slist
    import Samplelist2017_Nano14Dec2018 as Slist


IntLumiyear = {
    2016 : 35.92,
    2017 : 41.53,
    2018 : 59.74
}

YearInfo = {
    2016:{
        "IntLumi": 35.92,
        "GravitonMasses":[],
        "RadionMasses":[]
        },
    2017:{
        "IntLumi": 41.53,
        "GravitonMasses":[],
        "RadionMasses":[]
        },
    2018:{
        "IntLumi": 59.74,
        "GravitonMasses":[],
        "RadionMasses":[270, 320, 350, 400, 450, 500,  600, 650, 900]
        },
        
        }

IntLumi = YearInfo[Runyear]['IntLumi']
RadionMassPoints = YearInfo[Runyear]["RadionMasses"]


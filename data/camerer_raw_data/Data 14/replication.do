/* Replication of Ifcher and Zarghamee (2011) AER *//*                                                *//* Colin Camerer, Taisuke Imai, Gideon Nave       *//* 10/14/2015                                     *//*                                                */// Set memoryclearset mem 100m// Set working directorycd "SET WORKING DIRECTORY HERE"// Import dataclearuse all.dtareplace treatment = 0 if treatment == 2label define treatment 0 "Neutral" 1 "Positive"// Create variables as in IZg female=(gender==1)g arts_science=(college==1)g business=(college==2)g engineering=(college==3)g other=(college==4)g practicing=(practice_religion==1)g atheist_agnostic=(religion==1)g christian=(religion==3)g other_religion=(religion~=1&religion~=3)g asian=(race==5)g hispanic=(race==3)g white=(race==1)g other_race=~asian&~hispanic&~whiteg inc_lt100k=(family_income<=5)g inc_100_200k=(family_income==6|family_income==7)g inc_gt200k=(family_income==8|family_income==9)g time1=(delay==1)g time3=(delay==3)g time7=(delay==7)g time14=(delay==14)g time28=(delay==28)g time56=(delay==56)g m1134=(fv<15)g m1831=(fv>15&fv<20)g m2428=(fv>20&fv<30)g m3284=(fv>30&fv<40)g m5171=(fv>40)// Replication Table 3 Columns 1 and 5xi: reg pv treatment i.fv*i.delay, cluster(id)xi: reg pv treatment i.fv*i.delay if pv<fv, cluster(id)xi: reg pv treatment i.fv*i.delay i.college i.female i.race i.religion i.practicing i.family_income i.happiness if pv<fv, cluster(id)// Additional analysis//  video_happiness: "Did seeing the video clip make you" 1=Happier, 2=Neither happier, nor sadder, 3=Sadder//  video_mood: "Did seeing the video clip put you in" 1=A better mood, 2=Neither a better, nor a worse mood, 3=A worse moodtabulate video_happiness, gen(vh)tabulate video_mood, gen(vm)gen treat_vh1 = treatment*vh1gen treat_vh3 = treatment*vh3gen treat_vm1 = treatment*vm1gen treat_vm3 = treatment*vm3xi: reg pv treatment vh1 vh3 treat_vh1 treat_vh3 i.fv*i.delay if pv<fv, cluster(id)xi: reg pv treatment vm1 vm3 treat_vm1 treat_vm3 i.fv*i.delay if pv<fv, cluster(id)xi: reg pv treatment vh1 vh3 treat_vh1 treat_vh3 vm1 vm3 treat_vm1 treat_vm3 i.fv*i.delay if pv<fv, cluster(id)display _b[treatment] + _b[vh1] + _b[treat_vh1]test treatment + vh1 + treat_vh1 = 0display _b[treatment] + _b[vh3] + _b[treat_vh3]test treatment + vh3 + treat_vh3 = 0display _b[treatment] + _b[vm1] + _b[treat_vm1]test treatment + vm1 + treat_vm1 = 0display _b[treatment] + _b[vm3] + _b[treat_vm3]test treatment + vm3 + treat_vm3 = 0estimates table, varwidth(30) style(oneline)

xi: reg pv vh1 vh3 i.fv*i.delay i.college i.female i.race i.religion i.practicing i.family_income i.happiness if pv<fv & treatment == 0, cluster(id)xi: reg pv vm1 vm3 i.fv*i.delay i.college i.female i.race i.religion i.practicing i.family_income i.happiness if pv<fv & treatment == 0, cluster(id)xi: reg pv vh1 vh3 vm1 vm3 i.fv*i.delay i.college i.female i.race i.religion i.practicing i.family_income i.happiness if pv<fv & treatment == 0, cluster(id)xi: reg pv vh1 vh3 i.fv*i.delay i.college i.female i.race i.religion i.practicing i.family_income i.happiness if pv<fv & treatment == 1, cluster(id)xi: reg pv vm1 vm3 i.fv*i.delay i.college i.female i.race i.religion i.practicing i.family_income i.happiness if pv<fv & treatment == 1, cluster(id)xi: reg pv vh1 vh3 vm1 vm3 i.fv*i.delay i.college i.female i.race i.religion i.practicing i.family_income i.happiness if pv<fv & treatment == 1, cluster(id)xi: reg pv vh1 vh3 i.fv*i.delay i.college i.female i.race i.religion i.practicing i.family_income i.happiness if pv<fv, cluster(id)xi: reg pv vm1 vm3 i.fv*i.delay i.college i.female i.race i.religion i.practicing i.family_income i.happiness if pv<fv, cluster(id)xi: reg pv vh1 vh3 vm1 vm3 i.fv*i.delay i.college i.female i.race i.religion i.practicing i.family_income i.happiness if pv<fv, cluster(id)

// Manipulation checkkeep if question == 1bysort robin_died: tab robin_influencebysort treatment: tab video_happinessbysort treatment: tab video_moodregress video_happiness treatmentregress video_mood treatmentgen panas1 = (amusement+arousal+contentment+happiness+interest+relief+surprise)-(anger+confusion+contempt+disgust+embarrassment+fear+pain+sadness+tension)// Some subjects answered "0" but the result is robust if coding them as "1"'sreplace amusement = 1 if amusement == 0 replace arousal = 1 if arousal == 0replace contentment = 1 if contentment == 0replace happiness = 1 if happiness == 0replace interest = 1 if interest == 0replace relief = 1 if relief == 0replace surprise = 1 if surprise == 0replace anger = 1 if anger == 0replace confusion = 1 if confusion == 0replace contempt = 1 if contempt == 0replace disgust = 1 if disgust == 0replace embarrassment = 1 if embarrassment == 0replace fear = 1 if fear == 0replace pain = 1 if pain == 0replace sadness = 1 if sadness == 0replace tension = 1 if tension == 0gen panas2 = (amusement+arousal+contentment+happiness+interest+relief+surprise)-(anger+confusion+contempt+disgust+embarrassment+fear+pain+sadness+tension)ttest panas1, by(treatment)ttest panas2, by(treatment)
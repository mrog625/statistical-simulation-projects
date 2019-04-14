### Deck of Cards Shuffling Simulation Project

## 9-card deck exploration
m = 100; x = numeric(m) 
for(i in 1:m) 
{ 
  deck = sample(1:9, 9)     # Randomize deck of 9 cards 
  loc = match(1:9, deck)    # Index in 'deck' of each card 
  back = diff(loc) < 0      # Ts where rising seqs end (except last) 
  x[i] = sum(back) + 1      # Count rising sequences (including last) 
} 
mean(x); sd(x)
hist(x,breaks=(0:9+.5),col="wheat", xlab="Rising Sequences", 
     main="Simulated Distribution for a 9-Card Deck in Random Order") 


## Simulate 100000 random decks of cards and count the number of rising sequences for each
m = 100000; x = numeric(m) 
for(i in 1:m) 
  { 
    deck = sample(1:52, 52)   # Randomize deck of 52 cards 
    loc = match(1:52, deck)   # Index in 'deck' of each card 
    back = diff(loc) < 0      # Where rising seqs end (except last) 
    x[i] = sum(back) + 1      # Count rising sequences (including last) 
  } 
mean(x); sd(x)

# Percentage of random decks with 23 to 30 rising sequences
mean(x)+c(qnorm(.025),qnorm(.975))*sd(x)
mean(x)+c(qnorm(.005),qnorm(.995))*sd(x)
mean(x>=21 & x<=32)
mean(x>=23 & x<=30)

# Histogram of distirbution of rising sequences for random decks
cutp = seq(min(x)-1, max(x)) + .5 
hist(x, breaks=cutp, prob=T, col="wheat", xlab="Rising Sequences", 
     main="Simulated Distribution for a Deck in Random Order") 
curve(dnorm(x, 26.5, 2.1), add=T, col="blue", lwd=2) 


## Simulating riffle shuffles
par(mfrow=c(2,2)) 
shuffles=c(4,6,7,8)
num.mean.sd.prop=matrix(nrow=length(shuffles), ncol=4)
for (s in 1:length(shuffles))
{
m = 10000; x = numeric(m); n = shuffles[s]
for (i in 1:m) 
{ 
  deck = 1:52 # fresh deck 
  # Cut/Shuffle Deck n Times 
  for (j in 1:n) 
  { 
    cut.nr = rbinom(1, 50, 1/2) + 1   #Location of cut in deck
    ix = sort(sample(1:52, cut.nr))   #Randomly select where the cards-before-the-cut will be placed in the new deck
    nd = numeric(52)                  #Initialize temporary new deck variable
    nd[ix] = deck[1:cut.nr]           #Place the cards-before-the-cut in the new deck
    nd[nd==0] = deck[(cut.nr+1):52]   #Place the cards-after-the-cut in the new deck
    deck = nd                         #Redefine "deck" variable with the new deck's organization
  } 
  # Count Rising Sequences in Shuffled Deck 
  x[i] = sum(diff(match(1:52, deck)) < 0) + 1 
} 
mn = min(x); mx = max(x); cut = mn:(mx+1) - .5 
cut=10:40 - .5
hist(x, breaks=cut, col="wheat", prob=T, xlab="Rising Sequences",ylim=c(0,.4), 
     main=paste("Simulated", n, "Shuffles")) 
curve(dnorm(x, 26.5, 2.1), add=T, col="blue", lwd=2) 
num.mean.sd.prop[s,1]=shuffles[s];
num.mean.sd.prop[s,2]=mean(x); num.mean.sd.prop[s,3]=sd(x); 
num.mean.sd.prop[s,4]=mean(x>=21 & x<=32);
}
par(mfrow=c(1,1)) 
rbind(c("Shuffles","Mean Rising Seq","SD Rising Seq","% Random Decks"),num.mean.sd.prop)

# Calculate shared coverage
shared.coverage=numeric(52)
for (p in 1:52)
{
  shared.coverage[p]=min((sum(x==p)/length(x)), dnorm(p,26.5,2.1))
}
sum(shared.coverage)
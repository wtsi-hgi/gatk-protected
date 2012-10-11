pNonRef <- function(pls) {
  totalP = 10^(pls/-10)
  totalP = totalP / sum(totalP)
  pNonRef = sum(totalP[2:length(pls)])
}

toBi <- function(pls, second) {
  totalP = 10^(pls/-10)
  totalP = totalP / sum(totalP)
  
  if ( second ) {
    hom.x = c(1, 2, 3)
    het.alt = c(4,5)
    hom.alt = c(6)
  } else {
    hom.x = c(1, 4, 6)
    het.alt = c(2,5)
    hom.alt = c(3)
  }
  
  bi = c(sum(totalP[hom.x]), sum(totalP[het.alt]), sum(totalP[hom.alt]))
  bi = round(-10 * log10(bi))
  bi - min(bi)
}



print(toBi(c(0, 10, 20, 30, 40, 50), F))
print(toBi(c(0, 10, 20, 30, 40, 50), T))
print(toBi(c(0, 10, 10, 10, 10, 10), T))
print(toBi(c(0, 1, 2, 3, 4, 5), F))
print(toBi(c(0, 1, 2, 3, 4, 5), T))

print(toBi(c(0, 50, 50, 50, 50, 50), F))
print(toBi(c(0, 50, 50, 50, 50, 50), T))

print(toBi(c(50, 0, 50, 50, 50, 50), F))
print(toBi(c(50, 0, 50, 50, 50, 50), T))

print(toBi(c(50, 50, 0, 50, 50, 50), F))
print(toBi(c(50, 50, 0, 50, 50, 50), T))

print(toBi(c(50, 50, 50, 0, 50, 50), F))
print(toBi(c(50, 50, 50, 0, 50, 50), T))

print(toBi(c(50, 50, 50, 50, 0, 50), F))
print(toBi(c(50, 50, 50, 50, 0, 50), T))

print(toBi(c(50, 50, 50, 50, 50,  0), F))
print(toBi(c(50, 50, 50, 50, 50,  0), T))

print(toBi(c(5, 4, 3, 2, 1, 0), F))
print(toBi(c(5, 4, 3, 2, 1, 0), T))

print(toBi(c(10, 10, 10, 10, 10, 0), F))
print(toBi(c(10, 10, 10, 10, 10, 0), T))

# kep2cart
Algorithm implemented in Matlab R2020b that calculates Cartesian state vectors from a set of Keplerian orbit elements. 
The input to the algorithm is a set of Keplerian orbit elements: 
  The semi-major axis \texttt{a}, 
  eccentricity \texttt{e},
  argument of periapsis \texttt{w}, 
  longitude of ascending node \texttt{omega}, 
  inclination \texttt{i} 
  and the mean anomaly \texttt{M0}. 
Additionally the time \texttt{t} has to be specified. 
The outputs of the algorithm are the Cartesian position and velocity vectors \texttt{r} and \texttt{v}.

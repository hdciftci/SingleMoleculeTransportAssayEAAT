function FRET = efflux_function(time,opt_all)
monomer=3;
S=opt_all(1);
PEB1a=(1/6.02e23/(4/3*pi*opt_all(7)^3*10^(-24))*1e6);
V=(opt_all(8)*(opt_all(1)-0))./(opt_all(1)-0+opt_all(6)*(1+PEB1a/opt_all(4)))*monomer;
FRET=opt_all(1)/(opt_all(4)+opt_all(1))*opt_all(2)+opt_all(5);
for t=2:length(time)
    S=(S-V*opt_all(3));
    V=(opt_all(8)*S)./(S+opt_all(6)*(1+PEB1a/opt_all(4)))*monomer;
    FRET=[FRET;S/(S+opt_all(4))*opt_all(2)+opt_all(5)];
end
end






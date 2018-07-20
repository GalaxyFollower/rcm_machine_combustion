clear all
%Przypianie mechanizmu GRI-3.0, który zawiera m.in. złożone rekacje prowadzące do samozapłonu w gazach naturalnych. 
%Mechanizm GRI-3.0 był walidowany i optymalizowany także dla etanu. Parametry w reakcjach były optymalizowane,
%by otrzymać m.in. czas opóźnienia samozapłonu zgodny z badaniami doświadczalnymi.
gas=Solution('gri30.xml')

%Parametry powietrza  znajduącego się na zewnątrz komory spalania (parametry otoczenia)
T_ambient = 300.  % K
p_ambient = 1e5  % Pa
comp_ambient = 'O2:1, N2:3.76'

%Utworzenie rezerwuaru powietrza na zewnątrz komory spalania (otoczenia)
set(gas,'T', T_ambient,'P',p_ambient,'X',comp_ambient)
ambient_air = Reservoir(gas)

%Ustawienie parametrów początkowych mieszanki znajdującej się w cylindrze
%Odpowiednio temperatury, ciśnienia i składu mieszanki
%Utworzenie ubogiej mieszanki 1/lambda=0.52
set(gas,'T', 340.33,'P',2.97*oneatm(),'X','C2H6:0.52, O2:3.5, N2:13.16')

%Utworzenie reaktora, który stanowi komorę spalania
r1 = IdealGasReactor(gas)

%Ustawienie początkowej objętości cylindra i przypisanie jej do reaktora
vol_start=0.000907335
setInitialVolume(r1,vol_start)

%Utworzenie ruchomej ścianki tłoka oddzielającej objętość cylindra od otoczenia
piston = Wall(ambient_air, r1)
A=0.084*0.084/4*pi %Ustawienie powierzchni tłoka
setArea(piston,A)
	
%Ustawienie parametrów gazu składającego się z atomów wodoru
set(gas,'T', 500.0,'P',oneatm(),'X','H:1.0')

%Utworzenie rezerwuaru z atomami wodoru
%które wstrzyknięte do komory spalania mają modelować zapłon iskrowy 
igniter = Reservoir(gas)

%Przepływ masowy atomów wodoru z rezerwuaru do objętości cylindra jest
%modelowany zależną od czasu funkcją Gaussa
amplitude=0.0009; %z zapłonem
%amplitude=0.000; %brak zapłonu
t0 = 0.0293 % czas zapłonu, wstrzyknięcia atomów wodoru do komory spalania
m3=MassFlowController(igniter,r1);
mdot=gaussian(amplitude, t0, 0.0003);
setFunction(m3, mdot)

%Utworzenie połączeń między rezerwuarami i reaktorem, przypisanie komory spalania do symulacji
sim123 = ReactorNet([r1])

%Utworzenie pustych macierzy potrzebnych podczas obliczeń i do archiwizacji wyników 
Press=[]; %ciśnienie
Vol=[]; %objętość
Temp=[]; %temperatura
Tim=[]; %czas
dQdt=[]; %heat release rate, szybkość wydzielania ciepła
dQdt4=[]; %alternatywny sposób liczenia szybkości wydzielania ciepła
dpdt=[]; %różniczka ciśnienia po czasie
dVdt=[]; %różniczka objętości po czasie
Qdt=[]; % wydzielone ciepło
x_ub=[]; %udział masowy frakcji niespalonej
gamma=[]; %wykładnik adiabaty

massFuel0=massFraction(r1,'C2H6'); %masa etanu w cylindrze w chwili początkowej

%wartości prędkości tłoka
piston_velocity=[552.924257e-3 1207.673032e-3 2604.261203e-3 4096.098903e-3 6705.682123e-3 9753.71566e-3 10963.79976e-3 9161.470345e-3 11251.03518e-3 6276.102e-3 6096.5825e-3 3698.214118e-3 1395.119063e-3 -1522.798929e-3 -3658.656129e-3 -6136.002245e-3 -5806.578477e-3 -9580.858509e-3 -9792.539e-3 -9712.99375e-3 -8022.472364e-3];

%wartości czasów końcowych, do których trwa symulacja z zadną prędkością tłoka
t_end=[0.01043 0.01169 0.01381 0.01743 0.02132 0.02668 0.02874 0.02961 0.03071 0.03141 0.03201 0.03218 0.0325 0.03278 0.03309 0.03358 0.03509 0.0367 0.0397 0.03986 0.04277];

i=0; %pomocniczy iterator

t = 0.0014 %czas początkowy
setInitialTime(sim123, t); %przypisanie początkowego czasu
setVelocity(piston,polynom(piston_velocity(1))); %przypisanie prędkości tłokowi

%Proces symulacji wraz z obliczaniem oprócz podstawowych parametrów jak T,p,V,
%także HRR,wydzielone ciepło i udział frakcji spalonej.
%W jednej pętli while prowadzone są obliczenia tylko dla jednej prędkości tłoka
%i dla pewnego zakresu czasu.

while t < t_end(1)
	t_pocz=t;
	t=t+1.0e-6; %maksymalny krok czasowy 1e-6, solver może go zmniejszać w razie konieczności
	i=i+1;
	advance(sim123,t) %solver
	
	%Przypisanie wyników
	Press(i)=pressure(r1);
	Vol(i)=volume(r1);
	Temp(i)=temperature(r1);
	Tim(i)=t;

	if i==1 %przypisanie zerowych wartości w pierwszym kroku
		dQdt(i)=0;
		Qdt(i)=0;
		dQdt4(i)=0;

	else %obliczenia w dalszych krokach czasowych
		
		%obliczenia wykładnika adiabaty
		gamma(i)=cp_mole(gas)/cv_mole(gas); 
		
		%obliczenia różniczki objętości po czasie
		dVdt(i)=(Vol(i)-Vol(i-1))/(Tim(i)-Tim(i-1)); 
		
		%obliczenia różniczki ciśnienia po czasie	
		dpdt(i)=(Press(i)-Press(i-1))/(Tim(i)-Tim(i-1)); 
		
		%obliczenia szybkości wydzielania ciepła
		dQdt(i)=gamma(i)/(gamma(i)-1)*Press(i)*dVdt(i)+(1/(gamma(i)-1))*Vol(i)*dpdt(i); 
		
		%całkowanie HRR, by otrzymać wydzielone ciepło
		Qdt(i)=(dQdt(i)+dQdt(i-1))*0.5*(Tim(i)-Tim(i-1))+Qdt(i-1);
		
		%alternatywny sposób obliczenia HRR		
		dQdt4(i)=-Vol(i)*(cp_mass(gas)-cv_mass(gas))*mass(r1)*Temp(i)*enthalpies_RT(gas)'*netProdRates(gas); 
		end
	x_un(i)=(massFraction(r1,'C2H6'))/massFuel0; %obliczenie udziału masowego frakcji niespalonej
	
	end

%Analogicznie druga pętla
setInitialTime(sim123, t); %ustawienie końcowego czasu z poprzedniej pętli jako początkowego dla nowej
setVelocity(piston,polynom(piston_velocity(2))); %przypisanie drugiej prędkości tłoka
while t < t_end(2) %zmiana końcowego czasu na drugi
	t=t+1.0e-6;
	i=i+1;
	advance(sim123,t)
	
	Press(i)=pressure(r1);
	Vol(i)=volume(r1);
	Temp(i)=temperature(r1);
	Tim(i)=t;

	gamma(i)=cp_mole(gas)/cv_mole(gas);
	dVdt(i)=(Vol(i)-Vol(i-1))/(Tim(i)-Tim(i-1));
	dpdt(i)=(Press(i)-Press(i-1))/(Tim(i)-Tim(i-1));	
	dQdt(i)=gamma(i)/(gamma(i)-1)*Press(i)*dVdt(i)+(1/(gamma(i)-1))*Vol(i)*dpdt(i);
	Qdt(i)=(dQdt(i)+dQdt(i-1))*0.5*(Tim(i)-Tim(i-1))+Qdt(i-1);
	dQdt4(i)=-Vol(i)*(cp_mass(gas)-cv_mass(gas))*mass(r1)*Temp(i)*enthalpies_RT(gas)'*netProdRates(gas);
	x_un(i)=(massFraction(r1,'C2H6'))/massFuel0;
	end	
	
%Analogicznie utworzono pozostałe pętle dla t_end(3)...t_end(21) i prędkości od piston_velocity(3) do piston_velocity(21)	

%Na samym końcu wydrukowano minimalne i maksymalne wartości
minimumVol=min(Vol)/3 %1/3 minimalnej objętości [m^3]
maximumVol=max(Vol)/3 %1/3 maksymalnej objętości [m^3]
eps=max(Vol)/min(Vol) %stopień sprężania
maximumPress=max(Press)/10^6 %maksymalne ciśnienie [MPa]
maximumTemp=max(Temp) %maksymalna temperatura
Qmax=max(Qdt)/3 %1/3 maksymalnego wydzielonego ciepła [J]
HRRmax=max(dQdt)/3000 %1/3 maksymalej szybkości wydzielania ciepła [kW]

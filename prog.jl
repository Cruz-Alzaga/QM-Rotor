###################################################################
### GENERACIÓN DE CADENAS Y CÁCULO DE ACCIÓN Y CARGA TOPOLÓGICA ###
###################################################################

function action(phi,mi,t)
      s = 0.0

      for i = 1:t
            s += (phimod(phi[mod(i,t)+1]-phi[i]))^2.0
      end

      return (mi/2.0)*s
end


function tcharg(phi,t)
      q = 0.0

      for i = 1:t
            q += phimod(phi[mod(i,t)+1]-phi[i])
      end 

      return q/(2.0*pi)     
end


function phimod(α)
      α = mod(α,2.0*pi)
      if α > pi
            return α - 2.0*pi
      elseif α < -pi
            return α + 2.0*pi
      end

      return α
end


function update!(phi,mi,t)
      a = 0.0

      for i = 1:t # Bucle para los elementos de la configuración j
            ϵ = 2.0*(rand(Float64)-0.5) # Valor aleatorio en [-1,+1]
            phip = phi[i]+ϵ*3.0/sqrt(mi) # Elegimos el intervalo 3σ y normalizamos a 2π
            pacep = min(1.0,exp(-(mi/2.0)*((phimod(phi[mod(i,t)+1]-phip))^2+(phimod(phip-phi[mod(i-2+t,t)+1]))^2))/exp(-(mi/2.0)*((phimod(phi[mod(i,t)+1]-phi[i]))^2+(phimod(phi[i]-phi[mod(i-2+t,t)+1]))^2)))
            r = rand(Float64)
            # Aplicamos la selección de M-H
            if r < pacep
                  phi[i] = phip 
                  a += 1.0
            else 
                  a += 0.0
            end

            
      end

      # ahora un step de winding
      #a += windstep!(phi,mi,t)

      return a/(t )#+ 1.0)
end


function winding!(phi, t)
      for i = 1:t
          phi[i] = mod(phi[i] + (i-1)*2.0*pi/t , 2.0*pi)
      end
end
  

function antiwinding!(phi, t)
      for i = 1:t
          phi[i] = mod(phi[i] - (i-1)*2.0*pi/t , 2.0*pi)
      end
end


function windstep!(phi, mi, t)

      s_inicial = action(phi, mi, t)
      d = 0

      # Elegimos de manera aleatoria si hacemos wind o antiwind
      if rand() < 0.5
            winding!(phi, t)
            d += 1
      else
      antiwinding!(phi, t)
      end

      s_final = action(phi, mi, t)

      # Accept-reject de metropolis. Si se rechaza, deshacer el winding/antiwinding.
      pacep = exp(-(s_final - s_inicial))
      r = rand(Float64)

      # Aplicamos la selección de M-H
      if r > pacep
            if d==1
                  antiwinding!(phi, t)
            else   winding!(phi, t)                     
            end
            return 0.0
      else return 1.0
      end
      
end


function sq(mi,t,n,dis)
      phi = zeros(t)
      s = zeros(n)
      q = zeros(n)
      acep = 0.0

      #Damos un valores aleatorios a la primera configuración
      for i = 1:t
            phi[i] = 2*pi*rand(Float64)
      end 

      for i = 1:1000 # Bucle para alcanzar el equilibrio
            # Actualizamos la coordenada
            update!(phi,mi,t)
      end

      # Calculamos la acción y la carga topológica de la primera configuración
      s[1] = action(phi,mi,t)
      q[1] = tcharg(phi,t)

      for j = 2:n # Bucle de las distintas configuraciones
            # Actualizamos la coordenada
            acep += update!(phi,mi,t)
            # Calculamos la acción de la configuración final j
            s[j] = action(phi,mi,t)
            q[j] = tcharg(phi,t)

            if dis == true      
                  if rem(10*j,n) == 0 
                     println(j," ", s[j] ," ", q[j])
                  end 
            end
      end # Fin bucle de las distintas configuraciones

      acep = acep*100/n
      if dis == true
            println(raw"Acceptance rate: ",acep)
      end

      

      return s,q

end



###################################################################
######## CÁLCULO FUNCIÓN A DOS PUNTOS Y MASA EFECTÍVA #############
###################################################################

function twopf(mi,t,b,n)
      # Inicializamos a cero los vectores
      f = round(Int,n/b)
      calph = zeros(t,f)
      c = zeros(t)
      varc = zeros(t)
      dc = zeros(t)
      d = 0.0
      m = zeros(t)
      dm = zeros(t)
      phi = zeros(t)


      for i = 1:t
            phi[i] = 2*pi*rand(Float64)
      end 
      
      println("Inicializamos el calculo de la función a dos puntos")
      for α = 1:f
            for i = 1:b # Hacemos la media con b configuraciones
                  update!(phi,mi,t)
                  # Con una misma configuración calculamos cada t y empleamos la simetría de traslación
                  for j = 1:t # Calculamos c_{α} par acada t  
                        for k = 1:t # Aplicamos la simetría de traslación
                              calph[j,α] += phimod(phi[mod1(j+k-1,t)])*phimod(phi[k])                         
                        end
                  end
            end
            if mod(10*α,f) == 0
                  println("Completado: ", 100*α/f, "%")
            end
      end

      calph .= calph./(b*t)

      for j=1:t 
            for i = 1:f
                  c[j] += calph[j,i]
            end
            c[j] = c[j]/f

            for i = 1:f
                  varc[j] += (calph[j,i]-c[j])^2.0
            end

      end



      for i = 1:t
            phi[i] = 2*pi*rand(Float64)
      end 
      
      println("Inicializamos el calculo de la contribución de vacío")
      for i = 1:n
            update!(phi,mi,t)
            for j = 1:t
                  d += phimod(phi[j])
            end
            if mod(10*i,n) == 0
                  println("Completado: ", 100*i/n, "%")
            end
      end

      d = d/(n*t)
      c .= c .- (d^2.0)
      for i = 1:f
            dc[i] = sqrt(varc[i]/(f^2.0)+(4.0*d^2.0)/(n*t))
      end
      

      println("Inicializamos el calculo de la masa efectiva") 
      for i = 1:t-1
            m[i] = -log(abs(c[i+1]/c[i]))
            dm[i] = dc[i]*sqrt((1/c[i+1]^2+1/c[i]^2))
      end
      
      return c, dc, m, dm
end


function varc(mi,t,b,n)
      # Inicializamos a cero los vectores
      f = round(Int,n/b)
      calph = zeros(t,f)
      c = zeros(t)
      varc = zeros(t)
      phi = zeros(t)

      for i = 1:t
            phi[i] = 2*pi*rand(Float64)
      end 
      
      println("Inicializamos el calculo de la función a dos puntos")
      for α = 1:f
            for i = 1:b # Hacemos la media con b configuraciones
                  update!(phi,mi,t)
                  # Con una misma configuración calculamos cada t y empleamos la simetría de traslación
                  for j = 1:t # Calculamos c_{α} par acada t  
                        for k = 1:t # Aplicamos la simetría de traslación
                              calph[j,α] += phimod(phi[mod1(j+k-1,t)])*phimod(phi[k])                         
                        end
                  end
            end
            if mod(10*α,f) == 0
                  println("Completado: ", 100*α/f, "%")
            end
      end

      calph .= calph./(b*t)

      for j=1:t 
            for i = 1:f
                  c[j] += calph[j,i]
            end
            c[j] = c[j]/f

            for i = 1:f
                  varc[j] += (calph[j,i]-c[j])^2.0
            end

      end

      return varc
end

#plot([c,m],label=["C(t)" "M(t)"])
#=cm=scatter(c,label="Two-point function",yerror=dc,markershape= :cross, mc= :red)
plot!(m,label="Efective mass",yerror=dm,markershape= :xcross, mc= :green)=#



###################################################################
############### CALCULO DEL TIEMPO DE CORRELACIÓN #################
###################################################################

function bin(dat,b)

      n = length(dat)
      nb = round(Int,n/b)
      adat = 0.0
      s1 = 0.0
      σ = 0.0

      for i = 1:n
            adat += dat[i]
      end
      adat = adat/n

      for j = 1:n
            s1 += (dat[j]-adat)^2.0                        
      end
      s1 = s1/(n^2.0)

      for j = 1:nb
            θ = 0.0
            for k = 1:b
                  θ += dat[k+(j-1)*b]
            end
            θ = θ/b
            σ += (θ-adat)^2.0                        
      end
      σ = σ/(s1*nb^2.0)
      
      return σ

end


function binchain(dat,b)
      lb = length(b)
      σ = zeros(lb)

      for i = 1:lb
            σ[i] = bin(dat,b[i])
      end

      return σ      
end


function binning(dat)

      # Cálculo del binsize
      #=n = length(dat)
      b = Vector{Float64}()
      for i = 20:10:n/1000
            if mod(n,i) == 0 
                  push!(b,i)
            end
      end=#
      b = [20,25,30,40,50,60,75,100,125,150,200,250,300,375,400,500,600,750,1000,1200,1500,2000,2500,3000] #binning for N=300000  
      #b = [200,250,300,400,500,750,1000,1500,2000,3000,5000,6000,10000,15000,30000] #binning for N=3000000 19 
      #b = [10,20,25,40,50,100,125,200,250,500,1000] #binning for N=10000 11
      l = length(b)
      nb = zeros(l)
      nb .= 1.0./b
      
      σ = zeros(l)
      σ = binchain(dat,b)
      
      return nb, σ

end

#=
using DelimitedFiles
touch("bindat3.txt")
bindat = open("bindat3.txt","w")
writedlm(bindat,q)
close(bindat)
q = readdlm("bindat2.txt", '\t', Float64, '\n')

using Plots
bin=plot(ylabel="σ",xlabel="1/b")
mi=[1.0,2.0,3.0,4.0,5.0];
mi=[1.0,1.5,2.0,2.5,3.0];
mi=[0.5,1.0,1.5,2.0,2.5];
sig = zeros(5)
lis=[:solid, :dash, :dot, :dashdot, :dashdotdot];
for i = 1:5
      println("Iniciando el cáclculo para Î = "*string(mi[i]))
      b,s=binning(mi[i],round(Int,1000*mi[i]),300000)
      sig[i] = s[24]/2
      plot!(b,s,label=("Î = "*string(mi[i])),linestyles=lis[i])
end

binp
savefig(bin,"")

mi .= 1.0./mi.^2
pl=plot(mi,sig,ylabel="τ",xlabel="1/Î^2",label=false)
scatter!(mi,sig,label="Correlation time dependence on Î")
savefig(pl,"")
=#



###################################################################
#### COMPROBACIÓN DEL USO DE LA DISCRETIZACIÓN QUANTUM PERFECT ####
###################################################################

function chi(mi,t,n)
      _,q = sq(mi,t,n,false)
      exq2 = 0.0
      s = 0.0
      for i = 1:n
            exq2 += q[i]^2.0
      end

      exq2 = exq2/n

      for i = 1:n
            s = (q[i]^2.0-exq2)^2.0 
      end

      s = sqrt(s/n)

      return 2.0*exq2*mi/t, 2.0*s*mi/t
            
      
end


function chicp(mi)

      int = 0.0
      intn = 0.0 

      int, _ = quadgk(x -> (x^2)*exp(-mi*(x^2)/2.0), -pi, pi)
      intn, _ = quadgk(x -> exp(-mi*(x^2)/2.0), -pi, pi)
            
      return mi*int/(intn*2.0*pi^2)
end


function chis(mi)

      int = 0.0
      intn = 0.0 

      int, _ = quadgk(x -> (x^2)*exp(-mi*(1-cos(x))), -pi, pi)
      intn, _ = quadgk(x -> exp(-mi*(1-cos(x))), -pi, pi)
            
      return mi*int/(intn*2.0*pi^2)
end


function grachi(p)

      x = zeros(p)
      #int = zeros(p)
      #σ = zeros(p)
      intc = zeros(p)
      ints = zeros(p)
      exact = zeros(p)
      
      mii = 0.0
      mif = 5.0
      stp = (mif-mii)/(p-1)

      for i = 1:p
            
            x[i] = (i-1)*stp+mii
            #int[i], σ[i]  = chi((i-1)*stp+mii,100,n)
            intc[i] = chicp((i-1)*stp+mii)
            ints[i] = chis((i-1)*stp+mii)
            exact[i] = 1/(2.0*pi^2)

            if mod(10*i,p) == 0 
                  println(x[i]," ", intc[i], " ", ints[i])
            end
      end

      return x, intc, ints, exact
end 

#= using QuadGK

n = 20
xi = zeros(n)
int = zeros(n)
s = zeros(n)
mii = 0.0
mif = 5.0
stp = (mif-mii)/(n-1)

for i = 1:n
      xi[i] = (i-1)*stp+mii
      int[i], s[i]  = chi((i-1)*stp+mii,100,1000000)
      
      if mod(10*i,n) == 0 
            println(xi[i]," ", int[i])            
      end
end
      
pl=plot(x,[exact, intc, ints],label=["Exact solution" "Classically perfect" "Standard" ], linestyles=[:solid :dash :dashdot])
plot!(xlabel="ξ/2",ylabel="χ·ξ",legend=:bottomright)
scatter!(xi, int, yerror=s, label="Computational", markershape= :xcross, ms=4)
plot!(yrange=[0.0,0.1])
pl=plot(x,[exact, intc],label=["Exact solution" "Classically perfect"], linestyles=[:dash :dashdot])
=#

    
  
###################################################################
### GENERACIÓN DE CADENAS Y CÁCULO DE ACCIÓN Y CARGA TOPOLÓGICA ###
###################################################################

function action(phi,mi,t)

      s = 0.0
      for i = 1:t
            s += (phimod!(phi[mod(i,t)+1]-phi[i]))^2.0
      end

      return mi*s/2.0
end


function actionobc(phi,mi,t)

      s = 0.0
      for i = 1:t-1
            s += (phimod!(phi[i+1]-phi[i]))^2.0
      end

      return mi*s/2.0
end


function tcharg(phi,t)
      q = 0.0

      for i = 1:t
            q += phimod!(phi[mod(i,t)+1]-phi[i])
      end 

      return q/(2.0*pi)     
end


function phimod!(α)
      α = mod(α,2.0*pi)
      if α > pi
            return α - 2.0*pi
      elseif α < -pi
            return α + 2.0*pi
      end

      return α
end


function update!(phi,mi,t,acep)

      a = 0.0
      aw = 0.0
      s = 0.0
      sp = 0.0

      for i = 1:t # Bucle para los elementos de la configuración j
            ϵ = 2.0*(rand(Float64)-0.5)/4.0 # Valor aleatorio en [-1,+1]
            phip = phi[i] + ϵ*3.0/sqrt(mi) # Elegimos el intervalo 3σ y normalizamos a 2π
            s = (mi/2.0)*((phimod!(phi[mod(i,t)+1]-phi[i]))^2+(phimod!(phi[i]-phi[mod(i-2+t,t)+1]))^2)
            sp = (mi/2.0)*((phimod!(phi[mod(i,t)+1]-phip))^2+(phimod!(phip-phi[mod(i-2+t,t)+1]))^2)
            pacep = exp(-sp)/exp(-s)
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
      aw += windstep!(phi,mi,t)
      if acep == true
            return a/t, aw
      end

      return
end


function updateobc!(phi,mi,t,acep)

      a = 0.0
      s = 0.0
      sp = 0.0


      # Actualizar el valor φ(1)
      ϵ = 2.0*(rand(Float64)-0.5)/4.0 # Valor aleatorio en [-1,+1]
      phip = phi[1] + ϵ*3.0/sqrt(mi) # Elegimos el intervalo 3σ y normalizamos a 2π
      s = (mi/2.0)*((phimod!(phi[2]-phi[1]))^2)
      sp = (mi/2.0)*((phimod!(phi[2]-phip))^2)
      pacep = exp(-sp)/exp(-s)
      r = rand(Float64)
      # Aplicamos la selección de M-H
      if r < pacep
            phi[1] = phip 
            a += 1.0
      else 
            a += 0.0
      end

      for i = 2:t-1 # Bucle para los elementos de la configuración j
            ϵ = 2.0*(rand(Float64)-0.5)/4.0 # Valor aleatorio en [-1,+1]
            phip = phi[i] + ϵ*3.0/sqrt(mi) # Elegimos el intervalo 3σ y normalizamos a 2π
            s = (mi/2.0)*((phimod!(phi[mod(i,t)+1]-phi[i]))^2+(phimod!(phi[i]-phi[mod(i-2+t,t)+1]))^2)
            sp = (mi/2.0)*((phimod!(phi[mod(i,t)+1]-phip))^2+(phimod!(phip-phi[mod(i-2+t,t)+1]))^2)
            pacep = exp(-sp)/exp(-s)
            r = rand(Float64)
            # Aplicamos la selección de M-H
            if r < pacep
                  phi[i] = phip 
                  a += 1.0
            else 
                  a += 0.0
            end
      end

      # Actualizar el valor φ(T)
      ϵ = 2.0*(rand(Float64)-0.5)/4.0 # Valor aleatorio en [-1,+1]
      phip = phi[t] + ϵ*3.0/sqrt(mi) # Elegimos el intervalo 3σ y normalizamos a 2π
      s = (mi/2.0)*((phimod!(phi[t]-phi[t-1]))^2)
      sp = (mi/2.0)*((phimod!(phip-phi[t-1]))^2)
      pacep = exp(-sp)/exp(-s)
      r = rand(Float64)
      # Aplicamos la selección de M-H
      if r < pacep
            phi[t] = phip 
            a += 1.0
      else 
            a += 0.0
      end

      if acep == true
            return a/t
      end

      return
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
      aw = 0.0

      #Damos un valores aleatorios a la primera configuración
      for i = 1:t
            phi[i] = 2*pi*rand(Float64)
      end 

      for i = 1:1000 # Bucle para alcanzar el equilibrio
            # Actualizamos la coordenada
            update!(phi,mi,t,false)
      end

      # Calculamos la acción y la carga topológica de la primera configuración
      s[1] = action(phi,mi,t)
      q[1] = tcharg(phi,t)

      for j = 2:n # Bucle de las distintas configuraciones
            # Actualizamos la coordenada
            a , b = update!(phi,mi,t,true)
            acep += a
            aw += b
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
      aw = aw*100/n
      if dis == true
            println(raw"Acceptance rate: ",acep,"%")
            println(raw"Winding cceptance: ",aw,"%")
      end

      

      return s,q

end


function aceprate(mi,t)

      phi = zeros(t)
      for i = 1:t
            phi[i] = 2*pi*rand(Float64)
      end 
      # Termalización
      for i = 1:5000
            update!(phi,mi,t,false)
      end

      acep = 0
      for i in 1:10000
            a,_=update!(phi,mi,t,true)
            acep += a
      end

      return acep/10000
end



###################################################################
############### CALCULO DEL TIEMPO DE CORRELACIÓN #################
###################################################################

function bin(dat,b)

      n = length(dat)
      nb = round(Int,n/b)
      adat = 0.0
      σ = 0.0

      for i = 1:n
            adat += dat[i]
      end
      adat = adat/n

      for j = 1:nb
            θ = 0.0
            for k = 1:b
                  θ += dat[k+(j-1)*b]
            end
            θ = θ/b
            σ += (θ-adat)^2.0                        
      end
      σ = σ/(nb^2.0-nb)
      
      return σ

end


function binnorm(dat,b)
      lb = length(b)
      σ = zeros(lb)

      for i = 1:lb
            σ[i] = bin(dat,b[i])
      end
      σ .= σ./bin(dat,1)

      return σ      
end


# using LinearAlgebra
function linefit(x::Vector, y::Vector)
    n = length(x)
    ny = length(y)

    if n != ny
      return ERROR
    end

    A = [x ones(n)]
    b = y
    coeffs = A \ b  # Solve the linear system
    
    slope = coeffs[1]
    intercept = coeffs[2]
    
    return slope, intercept
end


function corr_time(v)

      l = (length(v)-1)
      tmax = min(l, 10001)
      Gamma = autocor(v,collect(0:tmax))
      ϵ = 0.0
      for i in 1:tmax
            if Gamma[i] <= ϵ
                  return sum(autocor(v, collect(0:i-1))) + 1/2
            end
      end

      return 0.0
end



###################################################################
######## CÁLCULO FUNCIÓN A DOS PUNTOS Y MASA EFECTÍVA #############
###################################################################

#=
using DelimitedFiles, Statistics , StatsBase
=#

function tpfpbc(conf,t,dis)

      meanphi = zeros(t)
      for i in 1:t
            meanphi[i] = mean(conf[:,i])
      end

      Ct = similar(conf)
      Ct .= 0.0
      for j = 1:t  
            for k = 1:t 
                  it = mod1(j+k-1,t)
                  Ct[:,j] .+= (conf[:,it].-meanphi[it]).*(conf[:,k].-meanphi[k])/t                       
            end
            if mod(10*j,t)==0 && dis == true
                  print(" · ")
            end
      end

      if dis == true
            println(" ")
      end

      return Ct
end


function ctfpbc(conf,t)
      
      meanphi = zeros(t)
      for i in 1:t
            meanphi[i] = mean(conf[:,i])
      end

      Ct = tpfpbc(conf,t,true)

      Ctmean = zeros(t)
      for i in 1:t
            Ctmean[i] = mean(Ct[:,i])
      end

      Ctstd = zeros(t)
      l = length(Ct[:,5])
      for i in 1:t
            ct = corr_time(Ct[:,i])
            st = std(Ct[:,i])
            Ctstd[i] = sqrt(2.0*ct)*st/sqrt(l)
            if mod(10*i,t)==0
                  print(" · ")
            end
      end

      return Ctmean, Ctstd
end


function ctfpbcbin(mi,t,nb,b)

      phi = zeros(t)
      for i = 1:t
            phi[i] = 2*pi*rand(Float64)
      end 
      # Termalización
      for i = 1:40000
            update!(phi,mi,t,false)
      end

      conf = zeros(b,t)
      Ct = zeros(nb,t)
      for i in 1:nb
            for j in 1:b
                  update!(phi,mi,t,false)
                  conf[j,:] .= phi
            end
            
            Ctconf = similar(conf)
            Ctconf= tpfpbc(conf,t,false) 
            for j in 1:t
                  Ct[i,j] = mean(Ctconf[:,j])
            end

            if mod(10*i,nb)==0
                  print(" · ")
            end
      end
      println(" ")
      
      Ctmean = zeros(t)
      for i in 1:t
            Ctmean[i] = mean(Ct[:,i])
      end

      Ctstd = zeros(t)
      for i in 1:t
            Ctstd[i] = std(Ct[:,i])/sqrt(nb)
            if mod(10*i,t)==0
                  print(" · ")
            end
      end

      return Ctmean, Ctstd
end


function tpfobc(conf,t,t0,dis)

      meanphi = zeros(t)
      for i in 1:t
            meanphi[i] = mean(conf[:,i])
      end

      Ct = similar(conf)
      Ct .= 0.0
      for j = 1:t  
            Ct[:,j] .+= (conf[:,j].-meanphi[j]).*(conf[:,t0].-meanphi[t0])                      
            if mod(10*j,t)==0 && dis == true
                  print(" · ")
            end
      end

      if dis == true
            println(" ")
      end

      return Ct

end


function ctfobc(conf,t,t0)
      
      Ct = tpfobc(conf,t,t0,true)
      Ctmean = zeros(t)
      for i in 1:t
            Ctmean[i] = mean(Ct[:,i])
      end

      Ctstd = zeros(t)
      l = length(Ct[:,5])
      for i in 1:t
            ct = abs(corr_time(Ct[:,i]))
            st = std(Ct[:,i])
            Ctstd[i] = sqrt(2.0*ct)*st/sqrt(l)
            if mod(10*i,t)==0 
                  print(" · ")
            end
      end

      return Ctmean, Ctstd
end


function ctfobcbin(mi,t,t0,nb,b)

      phi = zeros(t)
      for i = 1:t
            phi[i] = 2*pi*rand(Float64)
      end 
      # Termalización
      for i = 1:40000
            updateobc!(phi,mi,t,false)
      end

      conf = zeros(b,t)
      Ct = zeros(nb,t)
      for i in 1:nb
            for j in 1:b
                  updateobc!(phi,mi,t,false)
                  conf[j,:] .= phi
            end
            
            Ctconf = similar(conf)
            Ctconf= tpfobc(conf,t,t0,false) 
            for j in 1:t
                  Ct[i,j] = mean(Ctconf[:,j])
            end

            if mod(10*i,nb)==0
                  print(" · ")
            end
      end
      println(" ")
      
      Ctmean = zeros(t-2*t0+1)
      for i in 1:t-2*t0+1
            Ctmean[i] = mean(Ct[:,i+t0-1])
      end

      Ctstd = zeros(t-2*t0+1)
      for i in 1:t-2*t0+1
            Ctstd[i] = std(Ct[:,i+t0-1])/sqrt(nb)
            if mod(10*i,t)==0
                  print(" · ")
            end
      end

      return Ctmean, Ctstd
end


function msfobc(Ct,Ctstd,d)

      Ms = zeros(d)
      Msstd = zeros(d)
      for i in 2:d
            Ms[i] = log(abs(Ct[i-1]/Ct[i+1]))/2.0
            Msstd[i] = sqrt((Ctstd[i+1]/Ct[i+1])^2+(Ctstd[i-1]/Ct[i-1])^2)/2.0
      end
      
      return Ms, Msstd
end


function msfpbc(Ct,Ctstd,d,t)

      Ms = zeros(d)
      Msstd = zeros(d)
      for i in 2:d
            f = abs((Ct[i+1]+Ct[i-1]+Ct[t-i+1]+Ct[t-i-1])/(Ct[i]+Ct[t-i]))
            Ms[i] = acosh(f)
            p = sqrt(1/(f^2-1))
            Msstd[i] = p*sqrt((Ctstd[i+1]/(2*Ct[i+1]))^2+(Ctstd[i-1]/(2*Ct[i-1]))^2+(Ctstd[t-i+1]/(2*Ct[t-i+1]))^2+(Ctstd[t-i-1]/(2*Ct[t-i-1]))^2+(Ctstd[i]*f/Ct[i])^2+(Ctstd[t-i]*f/Ct[t-i])^2)
      end
      
      return Ms, Msstd
end



###################################################################
#################### SUSCEPTIBILDAD TOPOLÓGICA ####################
###################################################################

function chi(mi,t,n)

      _,q = sq(mi,t,n,false)
      exq2 = 0.0
      s = 0.0
      
      q .= q.^2.0
      exq2 = mean(q)

      s = std(q)
      s = 2.0*corr_time(q)*s/sqrt(n)

      return 2.0*exq2*mi/t, 2.0*s*mi/t
end


function chiobc(mi,t,n)

      phi = zeros(t)
      for i = 1:t
            phi[i] = 2*pi*rand(Float64)
      end 
      println("Termalización") 
      for i = 1:40000
            updateobc!(phi,mi,t,false)
      end

      conf = zeros(n,t)
      println("Generamos las configuraciones")
      for i in 1:n
            updateobc!(phi,mi,t,false)
            conf[i,:] .= phi
            if mod(10*i,n)==0
                  print(" · ")
            end
      end
      println(" ")

      q = similar(conf)
      for i in 1:t
            tp = mod(i,t) + 1
            q[:,i] .= phimod!.(conf[:,tp]-conf[:,i])./(2*pi)
      end

      println("Cálculo de la función a dos puntos")
      #t0 = floor(Int,t/2)
      t0 = 5
      Ct = zeros(t)
      Ctstd = zeros(t)
      Ct , Ctstd = ctfobc(q,t,t0)

      return 2*Ct[t0]*mi, 2*Ctstd[t0]*mi
end


function plchi(n,p)
      
      xi = zeros(n)
      int = zeros(n)
      s = zeros(n)
      mii = 0.0
      mif = 5.0
      stp = (mif-mii)/(n-1)

      for i = 1:n
            xi[i] = (i-1)*stp+mii
            int[i], s[i]  = chi((i-1)*stp+mii,100,p)
      
            if mod(10*i,n) == 0 
                  println(xi[i]," ", int[i])            
            end
      end

      return xi, int, s
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
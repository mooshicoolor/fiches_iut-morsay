disp("Exercice 1.1:")
M=[0,1,1;0,0,1;1,0,0]
N=[0,0,1;1,1,1;1,0,1]
disp(modulo(M+N,2),"M+N=")
disp(modulo(M*N,2),"M.N=")


disp("Exercice2.2")
disp("f est injective si et seulement si pour chaque image sur (Z/2Z)7, il existe au plus un antécédent sur (Z/2Z)3 . Soit deux matrice M et M’ appartenant à (Z/2Z)3 ")
disp("Or deux matrices sont différentes lorsqu’il existe au moins un coefficient d’une matrice qui est différent de celui de l’autre matrice. De plus M et M’ sont différentes, donc leurs images sont différentes.")
disp("Ainsi pour deux antécédents différents de (Z/2Z)3 , on a deux images différentes de f. Donc f est injective.")
disp("Exercice2.3")
P=[[0, 1, 1];[ 1, 1, 0]; ones(2,3)]
G=[eye(3,3);P]
disp(P,"Matrice P:")
disp(G,"Matrice generatrice G:")
disp ("Ensemble des elements de C:")
function y=f(M)
    y=modulo(G*M,2)
endfunction
for x=0:1
for y=0:1
for z=0:1
disp(f([x;y;z]))
end
end
end

disp("Exercice 3.4")
H=[P,eye(4,4)]
disp(H,"H=")
disp(modulo(H*G,2),"H.G=")
disp("Exercice 3.5")
disp("Si c appartient à C alors c=G.[m1;m2;m3] donc H.c=H.G.[m1;m2;m3] Or H.G=0")


disp("Exercice 3.6")
disp("Si B recoit le message sans l`erreur alors il recevra r=c+e avec e=0")
disp(" donc il recevra c")
disp("Dans ce cas le produit de H par r est 0, sinon le message contient des erreurs")

disp("Exercice 4.7")
c1=[1;0;0;0;1;1;1]
c2=[1;0;1;1;1;0;0]
disp(modulo(c1+[1;0;1;1;0;1;1],2),"c1+[1;0;1;1;0;1;1]=")
disp(modulo(c2+[1;0;0;0;0;0;0],2),"c2+[1;0;0;0;0;0;0]=")
disp("On a bien une égalité")

disp("Exercice 4.8")
function y=S(M); y=modulo(H*M,2);endfunction
Y=S(c1)
disp(Y, "H.c=")
disp("S(c + e)=H.(c+e)=H.c+H.e=0+H.e")
disp("S(e)=H.e on a donc  bien S(c + e)=S(e)")
disp("Exercice 4.9")
disp("H appartient à M d"'indice 4,7 (4 ligne et 7 colonne) et M la matrice antécédente appartient à M d"'incide 7,1 ")
disp("(7 ligne et 1 colonne); Donc le produit des 2 matrices appartiendra à M d"'indice 4,1, il aura donc 4 coefficient appartenant à Z/2")
disp("Exercice 4.10")
function [colonne1, colonne2]=syndrome(H)
    nbligne=size(H,1)
    disp(nbligne,"nbligne de H=")
    nbelement=(2^nbligne)
    for a=(1:nbligne)
        rempli=0;
        for b=(1:2^a)
            colonne1(rempli+1:rempli+nbelement/(2^a),a)= modulo(b,2)
            rempli=((nbelement/(2^a))*b)
        end
    end
    nbcolonne=size(H,2)
    colonne2=zeros(nbelement,nbcolonne);
    rempli=1
    poids=1
    position=nbcolonne+1
while (rempli<nbelement & poids<nbcolonne)
    cpt=position-1
    while (rempli<nbelement & (cpt>0))
erreur=zeros(nbcolonne,1)
erreur(cpt,1)=1
        if(poids>1) then
            for (a=1:(poids-1))
            erreur(position-a+1,1)=1
            
        end
    end
            disp(erreur,"erreur=")
synd=modulo(H*erreur,2)
exist=%F
compteur3=1
while ((~exist) & (compteur3<nbelement))
    exist=(isequal(colonne1(compteur3,:),synd') & (isequal(colonne2(compteur3,:),zeros(1,nbcolonne))))
    compteur3=compteur3+1
end
if(exist) then
colonne2(compteur3-1,:)=erreur'
rempli=rempli+1
end
    cpt=cpt-1
end
if (position==2 | position==nbcolonne+1) then
        poids=poids+1
        position=nbcolonne
    else 
        position=position-1; 
end 
end 
endfunction
[X,Y]=syndrome(H)
disp(X,Y)
disp("Exercice 4.11")
disp("r=c+e d"'après l"'énoncé")
disp("r=c+e <=> r+e=c+e+e <=> r+e=c+2e")
disp("Or 2e est vecteur dont les coefficient sont de forme 2k et appartiennent à Z/2Z")
disp(" et 2k===0 [2] donc 2e est une matrice nulle")
disp("ainsi r+e=c")
disp("Exercice 4.12")
function M=correction(r)
    [SYNDROME,ERREUR]=syndrome(H)
    A=S(r);
    exist=%F
    cpt=1
    while (~exist & cpt<=16)
        exist=(isequal(SYNDROME(cpt,:),A'))
        cpt=cpt+1
    end
    if (exist) then
        e=ERREUR(cpt-1)
            c=modulo(r+e,2)
    else
        disp("introuvable")
        c=r
    end

    for a=(1:3)
        rempli=0;
        for b=(1:2^a)
            Mtemp(rempli+1:rempli+((2^3)/(2^a)),a)= modulo(b,2)
            rempli=((2^3/(2^a))*b)
        end
    end
    trouve=%F
    cpt=1
    disp(Mtemp(cpt,:),"1=")
    while (~trouve & cpt<2^3)
        F=f(Mtemp(cpt,:)')
        trouve=(isequal(c,F))
        cpt=cpt+1
    end
    if (trouve) then
        M=Mtemp(cpt-1,:)'
    else
        disp("non trouvé")
        M=zeros(3,1)
    end
endfunction
disp (correction([0;0;1;0;1;0;0]), "Le message est:")

disp("Exercice 4.15")

function G=transforme(G1)
    nbligne=size(G1,1)
    nbcolonne=size(G1,2)
    final=[eye(3,3);P]
    disp(G1,"G1=")
    disp(final,"Resultat attendu=")
    if (nbligne==7 & nbcolonne==3) then
    G=G1
    cpt=2
    while (~isequal(G(:,1),final(:,1)) & cpt<=3)
        G(:,1)=modulo(G(:,1)+G(:,cpt),2)
        cpt=cpt+1        
    end

    cpt=1
        while (~isequal(G(:,2),final(:,2)) & cpt<=3)
        G(:,2)=modulo(G(:,2)+G(:,cpt),2)
        cpt=cpt+2    
    end
    cpt=1
        while (~isequal(G(:,3),final(:,3)) & cpt<3)
        G(:,3)=modulo(G(:,3)+G(:,cpt),2)
        cpt=cpt+1        
    end
    if (~isequal(G,final)) then
        disp("conversion impossible")
    end
   else
       disp("Erreur taille matrice")
    end
    
    
endfunction
disp("Exercice 4.16")
G1=[1,1,1;0,1,1;1,0,1;1,1,0;1,0,0;0,0,1;0,0,1]
disp(transforme(G1),"Resultat de transforme:")

disp("nous pouvons donc en deduire que C et C"' sont égales")
disp (H*[0;0;1;0;1;0;0])

disp("fin du programme")

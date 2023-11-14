author : Kraemer Valentin 
date : 12 nov 2023
desc : description des problèmes abordés dans ce dossier


I. Problème d'advection diffusion stationnaire :

On se pose le problème 1D suivant : 

    Soient $\Omega = ]x_a, x_b[$ un ouvert borné de $\mathbb{R}$,   $\eta,~ b \geq 0,~\nu>0, ~ \bar{f}\in L^2(\Omega),~ u_a, ~ u_b \in \mathbb{R}$, on cherche à résoudre le problème suivant :
        trouver $\bar{u}\in H^2{\Omega}$ tel que :

        $$
            \eta u + b u' - \nu u'' = f sur \Omega 
            u(x_a) = u_a
            u(x_b) = u_b
        $$

On se propose de résoudre ce problème avec l'algorithme de Schwarz avec conditions de Robin à l'interface. On y ajoutera par la suite une surcouche de parallélisation avec la bibliothèque mpi4py de python ainsi qu'un terme en temps. 
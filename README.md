Author : Kraemer Valentin 
Date : 12 nov 2023 
Desc : description du dossier


# Organisation des dossiers 

Le github s'organise de la manière suivante : 

- un dossier params : les paramètres du problème 

- un dossier src : les sources (python, matlab), n'étant modifiés que dans la branche test ou après validation COMPLETE 

- un dossier trash : des brouillons en tous genres

- un dossier DEMO : les démos sous forme de notebooks ou de latex


# Organisation des branches 
    - main
    - test
    - paral

# Description du problème 

On se donne pour le moment le problème d'advection diffusion stationnaire 1D. 
Schéma cell-centered VF, décentré hybride avec DD Robin 

# Etat actuel du projet : 

Code en orienté objet du schéma stationnaire avec condition de Dirichlet et Robin pour une solution monodomaine

Validation graphique de Dirichlet -> phase de validation en ordre à venir

Pour la suite : 

 - validation Robin

 - utilisation de MPI pour la parallélisation DD en exploitant l'orienté objet pour éviter les problèmes d'indexation 


"""
    climato_tp(ds; fct::Function=sum, iduree=3, lead=1)

Calcul de la somme des précipitations pour tous les mois pour la durée spécifié par iduree. En ce moment, l'overlap des années n'est pas pris en charge: si iduree=3, on aura un Cube qui contient les accumulations sur 3 mois pour toutes les séquences de 3 mois, allant de janvier (janvier-février-mars) à décembre (décembre-janvier-février).
"""
function climato_tp(ds; fct::Function=sum, iduree=3, lead=1)
    # On veut éviter l'overlap

    moiscirc = CircularArray(1:12) 

    # timestamp = Array{DateTime}(undef, length(yearvec), length(mois))
    results = Array{YAXArray}(undef, length(moiscirc))
    
    i=1 # compteur pour les résultats mensuels
    for imois in 1:12

        month1 = moiscirc[imois+lead]
        month2 = moiscirc[imois+iduree-1+lead]

        # On subset les mois voulu        
        obs_subset = subsample(ds, month1 = month1, month2=month2)
        
        # on applique la fonction voulu sur le subset suivant selon le temps
        results[i] = yearly_clim(obs_subset, fct=fct)
        i += 1
    end
    
    # On concatenate les cubes de la liste results
    return concatenatecubes(results, Dim{:Mois}(1:12))
    
end

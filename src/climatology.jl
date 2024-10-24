"""
    climato_tp(ds; fct::Function=sum, iduree=3, lead=0)

Calcul de la somme des précipitations pour tous les mois pour la durée spécifié par iduree. En ce moment, l'overlap des années n'est pas pris en charge: si iduree=3, on aura un Cube qui contient les accumulations sur 3 mois pour toutes les séquences de 3 mois, allant de janvier (janvier-février-mars) à octobre (octobre-novembre-décembre). 
"""
function climato_tp(ds; fct::Function=sum, iduree=3, lead=0)
    # On veut éviter l'overlap
    # mois = 1:12-iduree+1 # ancienne version sans lead
    mois = 1:12-iduree+1-lead

    # timestamp = Array{DateTime}(undef, length(yearvec), length(mois))
    results = Array{YAXArray}(undef, length(mois))
    
    i=1 # compteur pour les résultats mensuels
    for imois in mois    

        # On subset les mois voulu     
        # obs_subset = subsample(ds, month1 = imois, month2=imois+iduree-1)   # ancienne version sans lead
        obs_subset = subsample(ds, month1 = imois+lead, month2=imois+iduree-1+lead)
        
        # on applique la fonction voulu sur le subset suivant selon le temps
        results[i] = yearly_clim(obs_subset, fct=sum)
        i += 1
    end
    
    # On concatenate les cubes de la liste results
    return concatenatecubes(results, RangeAxis("Mois", 1:length(mois)))
    
end

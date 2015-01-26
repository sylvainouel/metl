#!/usr/bin/tclsh

# define the lists of concepts you will be using
#
set problem {*cell2switch}
set mouvement {*move}
set voisinage {*neighborhood}
set gain {*gain}
set liste_tabu {*tabu_list}
set generateur {*init_sol_gen}
#set croisement {*uniform_xover<problem> *single_point_xover<problem> *two_point_xover<problem>}
set croisement {*uniform_xover<problem>}
#set accept_scheme {*metropolis<problem,mouvement> *threshold_accept<problem,mouvement>}
set accept_scheme {*metropolis<problem,mouvement>}

set cooling_scheme {*cooling_geometric_steps<voisinage>}
set mutation {*cell2switch_mutation}
set ls {*no_localsearch<problem>}
#set selection {*select_random<problem> *select_rank<problem> *select_tournament<problem> *select_roulette<problem>}
set selection {*select_random<problem>}
#set remplacement {*replace_worst<problem> *replace_worst_from_pop<problem> *replace_worst_parent<problem>}
set remplacement {*replace_worst<problem>}




#############################
# lists of algorithm templates. Directly taken from the framework

set algo_voisinage {
    *descent<problem,mouvement,voisinage,generateur> 
    *descent_gain<problem,mouvement,gain,generateur> 
    *descent_fm<problem,mouvement,voisinage,generateur> 
    *descent_ns_omp<problem,mouvement,voisinage,generateur> 
    *descent_ns_mpi<problem,mouvement,voisinage,generateur>
    *simulated_annealing<problem,mouvement,voisinage,generateur,accept_scheme,cooling_scheme> 
    *tabu<problem,mouvement,voisinage,liste_tabu,generateur>
    *tabu_gain<problem,mouvement,gain,liste_tabu,generateur> 
    *tabu_ns_omp<problem,mouvement,voisinage,liste_tabu,generateur> 
    *tabu_ns_mpi<problem,mouvement,voisinage,liste_tabu,generateur> 
}

set algo_evolution {
    *evolution<problem,generateur,croisement,mutation,ls,selection,remplacement>
}

set algo_sa {
    *simulated_annealing<problem,mouvement,voisinage,generateur,accept_scheme,cooling_scheme> 
}


set algo [concat $algo_voisinage $algo_evolution]
set coop {*omp_blackboard_coop *omp_ring_coop *omp_reduce_coop *mpi_blackboard_coop *mpi_ring_coop *mpi_reduce_coop}


########### this functions makes recurcives substitutions from the
########### format in argument and the previously defined template
########### lists
proc expand {format} {
    # split format
    set algo_list {}
    set av [split $format " <>,"]
    set expand_something 0

    foreach v $av {
	if {$v=="" || [regexp {^\*} $v]} continue
	set expand_something 1
	
	global $v

	if {![info exists $v]} {
	    puts stderr "no subtitution found for $v, skiping"
	    continue
	}

	eval set x $$v
	foreach p $x {
	    if {$v==$p} continue		

	    if {![regsub "(\[^\*_a-zA-Z0-9])$v" $format "\\1$p" the_expand]} {
		if {![regsub "^$v" $format "$p" the_expand]} {
		    puts stderr "substritution {$v -> $p in $format} failed"
		    return
		}
	    }
	    set algo_list [concat $algo_list [expand $the_expand]]
	}
	break
    }
    if {$expand_something==0} {
	while {[regsub -all {>>} $format {> >} format]} {}
	lappend algo_list $format
    }
    return $algo_list
}


#### remove elemets matching regular expression exp from a list
proc l_filter {l exp} {
    set r [list]
    foreach i $l {
	if {![regexp $exp $i]} {
	    lappend r $i
	}
    }
    return $r
}



set algo_base [expand algo_voisinage]
set old_gen $generateur
set generateur $algo_base
set algo_evo_gen [expand algo_evolution]
set generateur $old_gen
set old_ls $ls
set ls $algo_base
set algo_evo_ls [expand algo_evolution]

set algo_list [concat $algo_base $algo_evo_gen $algo_evo_ls]

set algo_for_coop [l_filter $algo_list "^.descent"]
#puts $algo_for_coop
#crap
set algo_coop [expand coop<algo_for_coop>]

set algo_list [concat $algo_list $algo_coop]
puts stderr [llength $algo_list]
### faire varier les concepts possédant des implémentation génériques 1 facteur à la fois.
set ls {*no_localsearch<problem>}

set old_croisement $croisement
set croisement {*uniform_xover<problem> *single_point_xover<problem> *two_point_xover<problem>}
set l_crossover [expand algo_evolution]
set croisement $old_croisement


set old_accept $accept_scheme
set accept_scheme {*metropolis<problem,mouvement> *threshold_accept<problem,mouvement>}
set l_accept [expand algo_sa]
set accept_scheme $old_accept

set old_selection $selection
set selection {*select_random<problem> *select_rank<problem> *select_tournament<problem> *select_roulette<problem>}
set l_selection [expand algo_evolution]
set selection $old_selection

set old_remplacement $remplacement
set remplacement {*replace_worst<problem> *replace_worst_from_pop<problem> *replace_worst_parent<problem>}
set l_remplacement [expand algo_evolution]
set remplacement $old_remplacement


set problem {*qap_prob}
set generator {*qap_gen}
set mutation {*mutation}

set croisement {*path_xover_swap<problem,mouvement> *path_xover_insert<problem>}
set l_crossover2 [expand algo_evolution]


set algo_list [concat $algo_list $l_crossover $l_accept $l_selection $l_remplacement $l_crossover2]



### output algo_list in a compilable form
set i 0
foreach x $algo_list {
    regsub -all {\*} $x {} output
    set name "algo[incr i]"
    set output2 "{ $output $name; test_generator(&$name); }"
    if {[regexp {mpi_[^_]*_coop.*ns_mpi} $output]} {
	set output2 "// $output2 // NOT WORKING"
    }

    puts $output2
}

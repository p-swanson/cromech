----------------------------
***C***oa***R***se-grain ***O***pen-***M***olecul***A***r-m***ECH***anics
----------------------------
```
 ██████╗██████╗  ██████╗ ███╗   ███╗███████╗ ██████╗██╗  ██╗
██╔════╝██╔══██╗██╔═══██╗████╗ ████║██╔════╝██╔════╝██║  ██║
██║     ██████╔╝██║   ██║██╔████╔██║█████╗  ██║     ███████║
██║     ██╔══██╗██║   ██║██║╚██╔╝██║██╔══╝  ██║     ██╔══██║
╚██████╗██║  ██║╚██████╔╝██║ ╚═╝ ██║███████╗╚██████╗██║  ██║
 ╚═════╝╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚═╝╚══════╝ ╚═════╝╚═╝  ╚═╝
```

                                    .+++====+**#*+:                                                                                                              
                                  .*:               +-                                                                                                           
     .%-..                        =.                  :*.                                                                                                        
     +.  .=+:.                   :-                     .*.                                                                                                      
     =        =+:.               =         =@             .=                                                                                                     
     -           .+=:.          .+                         ::                                                                                                    
     -               .==-.      :=                       .==                                                                                                     
     -                   .-=:   :=                      .+                                                                                                       
     -                      .*. .=                       .#                                                                                                      
     -                        .#-=                        :.                                                                                                     
     .*:.                       .#%:                     :*.                                                                                                     
           :%=..                   .#*..  :%%+:........:*.#..                                                                                                    
               .++..                    .=#:..              :#:.                                                                                                 
                   .-#*--.. ..                .#+-.   ....     .+=.                                                                                              
                           -%-=-.                  -++:   ==.     +.                                                                                             
                           +.+=: .==#+.             .=      +      -.                                                                                            
                         .*.-=       .+..===+:      :.      -:=++.  .=                                                                                           
                         .: .*        .-++    .:--**#.      +    ..=..=.                                                                                         
                         =.  .+           .-:    .==#       =      :.  .=                                                                                        
                        ::    .#.           :=    .+      -#.      +.   .+.                                                                                      
                        =      .=            -...=:.     .#. .....       -.                                                                                      
                       :=       --                     .*.               .#                                                                                      
                       :=        +                    .+.                 @                                                                                      
                       :=        .+=:            .::-*-                   @                                                                                      
                       .+            ...........                          @                                                                                      
                        =                                                 @                                                                                      
                        =.                                                @                                                                                      
                        ..-=              .:                             :@                                                                                      
                           -***#+...+*+-:...-          :-.=:  .=.=-     .=.                                                                                      
                                            .%=...-#***.   ....   ..%: .=                                                                                        
                                                                      .:                                                                                         
                                                                                                                                                                 
                                                                             
                                                                                                                                          
**figure 1. Cro-Magnon carrying a club. circa 2021 AD.**
----------------------------

This is a repo for coarse grained (CG) openMM code...

insipred (copied) by examples by Davit Potoyan.. who you should google.... (link coming)

Essentially, a library for doing coarse-grained simulations for intrinsically disordered proteins.

I will continue to update, but for now this repo is to hold my code that I made 2 years ago... gosh, has it been two years? I am such a loser..

----------------------------

for reasons that should be obvious to me, but remain elusive, if you run your CG simulations on GPU (with 0 net charge) there is an incredible speed up... especially in parrallel. 
CPU-heavy softwares (i.e. LAMMPS, which is really, really good..) even when run on soley GPU are unable to perform to this degree (2500 atoms, 1ns per MINUTE on 4 L40s lol.. (data incoming) ) . hence, I am using openMM, BUT. openMM was designed for atomistic simulations on proteins, NOT CG particles.
one of Straubs students (also look him up) made a blog post about using openMM for CG simultations on lipids. but using openMM for CG simulations requires like 72 lines of code.. (and LAMMPS: 83 lines of code.. again, its really good. and i would perfer using it)  so i am creating a library for code that can set up a simulation (generate toplogies with the option for user-input potentials) in 1 line of code.
---------------                                                                                                                       
> [!WARNING]
>  this code was made in my free-time and is not affiliated with anyone or aything besides myself
----------------------------
I wil create an easy-to-use library that, uh, allows you to do CG simulations using OpenMM, taking advantage of its ability to use custom potentials ... okay i am getting bored now.. I will come back to this , but uh.. why are you still reading
> [!WARNING]
>  I did not use any AI in creating this.. i think that true knowledge comes from suffering (something a machine is incapable of)


                                                            

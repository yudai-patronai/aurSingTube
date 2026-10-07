<!--mermaid-attrs:splitRatio=39-->
```mermaid

graph LR

    classDef molecule fill:#E8F5E9,stroke:#388E3C,stroke-width:6px,font-size:24px;
    classDef reaction fill:#FFEB3B,stroke:#FBC02D,stroke-width:6px,shape:rect,font-size:22px;

    %% Species
    A((A)):::molecule
    A_((A*)):::molecule
    AA_((AA*)):::molecule
    Pp((Ppase)):::molecule
    A_pp((A_ppase)):::molecule
    
  

    %% Reactions
    A --> Rxn1[kcat]:::reaction --> A_
    R2 --> Rxn1
    
    R3 --> Rxn2[]:::reaction --> R4
    R1 --> Rxn2
    
    R4 --> Rxn3[Rxn 3]:::reaction --> R5
    R4 --> Rxn4[Rxn 4]:::reaction --> R6
    
   
```
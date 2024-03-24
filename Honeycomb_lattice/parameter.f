
      parameter (nlatt=10, nvol=2*nlatt**2)
      
      integer campo
      real :: J_1, J_2, J_3
      
      common/reticolo/campo(nvol)

      ! Primi vicini
      common/vicinip/nn1(nvol),nn2(nvol),nn3(nvol)

      ! Secondi vicini
      common/vicinis/nnn1(nvol),nnn2(nvol),nnn3(nvol),
     &               nnn4(nvol),nnn5(nvol),nnn6(nvol)

      ! Terzi Vicini
      common/vicinit/nnnn1(nvol),nnnn2(nvol),nnnn3(nvol)

      common/vario/bext,J_1,J_2,J_3
      
C================================================================
C file che permette di cambiare più facilmente i parametri della
C simulazione facendolo una sola volta qui piuttosto che diverse
C volte all'interno del codice della simulazione; fare un ciclo
C sui reticoli allocando e deallocando il campo vorrebbe dire una
C simulazione molto lunga con il rischio di non poter verificare
C che tutto sia funzionando bene e quindi potenziale perdita di
C tempo
C================================================================

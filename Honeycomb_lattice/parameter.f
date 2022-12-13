
      parameter (nlatt=40, nvol=2*nlatt**2)
      
      integer campo
      !indicizzati da zero per comodità
      common/reticolo/campo(0:nvol)
      common/vicini/n1(0:nvol),n2(0:nvol),n3(0:nvol)
      
      common/vario/bext,J_acc
      
C================================================================
C file che permette di cambiare più facilmente i parametri della
C simulazione facendolo una sola volta qui piuttosto che diverse
C volte all'interno del codice della simulazione; fare un ciclo
C sui reticoli allocando e deallocando il campo vorrebbe dire una
C simulazione molto lunga con il rischio di non poter verificare
C che tutto sia funzionando bene e quindi potenziale perdita di
C tempo
C================================================================

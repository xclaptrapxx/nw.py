Algorytm Needleman-Wunscha do globalnego dopasowania sekwencji

Opis
Program implementuje algorytm Needleman-Wunscha w języku Python, umożliwiając globalne dopasowanie dwóch sekwencji, takich jak DNA lub białka. Program przyjmuje dwie sekwencje zapisane w formacie FASTA i wyświetla wynik dopasowania.
Algorytm korzysta z macierzy punktacji i macierzy śledzenia, aby znaleźć najlepsze dopasowanie między dwiema sekwencjami, uwzględniając dopasowania, niedopasowania oraz kary za luki.

Wymagania

    Python 3.x
    Biblioteki:
        Biopython - do obsługi plików FASTA (pip install biopython)
        NumPy - do operacji na macierzach (pip install numpy)

Sposób użycia

    Upewnij się, że w katalogu głównym znajduje się plik seqs.fa z dwoma sekwencjami w formacie FASTA.
    Uruchom skrypt w środowisku Python:
    
    bash:
    python nw.py

    Program wczyta sekwencje z pliku np. seqs.fa i wyświetli wynik dopasowania oraz wyrównane sekwencje.

Przykładowe dane wejściowe

Poniżej znajduje się przykładowa zawartość pliku seqs.fa:

  >seq1
  GATTACA
  >seq2
  GCATGCU

Po uruchomieniu programu:
  Wynik dopasowania: 0
  Alignment 1: G-ATTACA
  Alignment 2: GCA-TGCU

Uwagi:

    Program wykorzystuje domyślną punktację: +1 za dopasowanie, -1 za niedopasowanie oraz karę -1 za lukę.
    Przed uruchomieniem upewnij się, że plik seqs.fa zawiera dwie sekwencje w formacie FASTA.   

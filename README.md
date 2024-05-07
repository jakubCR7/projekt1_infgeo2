## Informatyka Geodezyjna II
## Projekt 1 - Transformacje 
Program umożliwia implementację transformacji geodezyjnych zgodnie z potrzebami użytkownika.

## Wymagania do obsługi programu:
Program został stworzony z użyciem programu Python 3.11.5, zaimporotwana została biblioteka NumPy, pozwalająca na wykonywanie obliczeń numerycznych i naukowych. Program został napisany w dla systemu operacyjnego Microsoft Windows 10 PRO i wyższych.

## Funkcje programu:
**Użytkownik wybiera spośród dostępnych elipsoid:**
'''
GRS80
WGRS84
Krasowski
'''
**Użytkownik wybiera spośród dostępnych transformacji:**
'''
 1 = X, Y, Z --> phi, lam, h 
 2 = phi, lam, h --> X, Y, Z 
 3 = X, Y, Z --> neu 
 4 = BL --> X2000, Y2000 
 5 = BL --> X92, Y92
'''
### Działanie programu:
Uruchamiając program Python wymagane jest również podanie nazwy pliku ze współrzędnymi wraz z jego rozszerzeniem (((tutaj trzeba o nagłówku i zerach wpisać))). Następnie, zgodnie z instrukcjami programu należy wybrać elipsoidę. Następnie, poprzez wpisanie liczby od 1 do 5, dokonujemy wyboru transformacji. Opcje, spośród których dokonujemy wyboru są podane powyżej w punkcje **Funkcje programu**.

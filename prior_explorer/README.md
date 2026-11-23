# Prior Explorer v3 — priory a fakes vazene podle truth

Samostatna nova verze GUI. Puvodni Machine.C, Machine3.C ani starsi GUI se neupravuji.

## Spusteni

Rozbalte slozku `prior_explorer_v3` do hlavniho adresare unfolding projektu a ponechte vsechny jeji soubory pohromade.

```bash
python3 -m pip install -r prior_explorer_v3/requirements.txt
python3 prior_explorer_v3/prior_explorer_v3.py
```

GUI potrebuje Tkinter (`python3-tk` na Ubuntu). Priprava pouziva vase funkcni prostredi ROOT + RooUnfold. Samotne GUI nepotrebuje PyROOT.

## Jedna nova priprava je nutna

Stare v1/v2 NPZ balicky obsahuji pouze reco spektrum fakes, nikoliv jejich vazbu na truth. Z takoveho spektra nelze zpetne zjistit, jak se zmeni po prevazeni truth pT.

1. Kliknete **Prepare from Machine.C**.
2. Vyberte puvodni nominalni `Unfolding/Machine.C` (pripadne nominalni Machine2.C), ze ktereho vznikl puvodni balicek. **Nepouzivejte Machine3.C s aplikovanymi trial vahami ani verzi s rucne upravenou FONLL funkci.** Doslo by k dvojimu prevazeni. Zname znacky `priorFactor`, `wPriorPt`, `PriorTrial.root` priprava odmitne; libovolnou rucni zmenu FONLL automaticky poznat nemuze.
3. Nastavte stejny D0 rozsah a poradi cut/reshape jako predtim. Vychozi `cutBefore=true` znamena cut, potom reshape.
4. Priprava v izolovane kopii makra zachyti kazde standardni volani `Fake()` a k jeho puvodni vaze ulozi skutecne `_mcJetPt`, pripadne `_mcD0Z` nebo `_mcLambda[iAng]`. Neprepocitava klasifikaci podle noveho odhadu.
5. Export overi, ze soucet kazde mapy pres truth kategorie presne reprodukuje puvodni `Hfakes()` v kazdem reco binu. Pokud struktura vaseho makra neodpovida znamym mistum plneni, priprava skonci s vysvetlenim misto ticheho pouziti nekompletni mapy.

Jedna priprava nacita vsechny tri centrality a osm modelu (jet pT, z, pet angularit a pTD), jako predchozi verze. Pro jiny D0 rozsah vytvorte dalsi balicek. ROOT se zastavi pred unfoldingem. Podrobny log je v `PriorExplorerBundles/prepare_.../prepare.log`.

Pouziva standardni projektovou strukturu: `config.h`, `Unfolding/Machine.C`, `Data/` a vstupy `Output/Output_ICS2.root`, `Output/Output_D02.root`, `Output/FONLL_D0Pt_0_31.root`, `Output/FONLL_D0Pt_0_10.root`. Binning i cuty prebere z makra/configu. Nestandardni relativni zavislosti mohou vyzadovat upravu pripravy.

Stare NPZ lze otevrit pro porovnani, ale pouze v rezimu **Fixed**. Aplikace to zretelne oznaci; chybejici truth mapy nedoplnuje odhadem.

## Fake weighting

| Volba | Chovani |
|---|---|
| Truth weights / keep outside | Vychozi. Fakes se znamym truth pT v rozsahu dostanou pT cast stejne vahy jako matched/misses. Mimo truth pT rozsah zustava vaha 1. |
| Truth weights / edge extension | Zname truth hodnoty mimo rozsah pouziji vahu nejblizsiho krajniho binu. Jde o explicitni predpoklad pokracovani vah, nikoliv extrapolaci FONLL nebo fitu. |
| Fixed | Puvodni v1/v2 model: fake sablona je pevna, pouze s volitelnym celkovym faktorem Fake normalization. |

Nezname/neplatne truth pT (nekladne nebo nefinite) zustava vzdy na vaze 1. U druhe promenne se za nedostupne povazuji nefinite hodnoty a sentinel <= -900.

### Dulezity detail pro 2D

Vaha se zapisuje jako `w(x,y) = C * g(x) * h(x,y)`, kde g meni jet pT a h druhou promennou vcetne korelace s pT.

- Je-li c=d=0, vaha zavisi jen na truth pT. Pouzije se i pro fakes, jejichz truth angularita/z je mimo rozsah nebo nedostupna.
- Pri nenulovych c,d a `keep outside` se mimo rozsah druhe osy ponecha `h=1`, ale dostupna pT cast `C*g` se stale pouzije. Mimo rozsah pT se ponecha cela vaha 1.
- Pri `edge extension` se kazda dostupna souradnice omezi na svuj krajni truth bin a pouzije se takto prodlouzena binova vaha. Pokud druha souradnice chybi, pouzije se pouze pT cast.
- Power-law vaha v pT sloupci bez nominalnich truth udalosti neni urcena (nulovy jmenovatel). Fakes v takovem sloupci se ponechaji na vaze 1; stejna konvence plati pro krajni prodlouzeni tohoto sloupce.

Tato konvence zajistuje plynuly prechod pri c,d -> 0. Rezimy s fakovymi vahami tedy meni fake tvar i normalizaci; nesnazi se umele zachovat procento fakes.

## Co zustava pevne

`A[i,j]` jsou vazene matched pocty. `T[j]` zahrnuje matched i misses. `R[i,j]=A[i,j]/T[j]` zustava pevne a zahrnuje ucinnost. Vahy jsou konstantni v truth binech.

Nova predikce je

```
reco = scale * (R @ trial_truth + fake_factor * weighted_fake_template)
```

`weighted_fake_template` se pocita z nove exportovanych map podle skutecnych truth souradnic fakes. Vychozi truth vaha se normalizuje tak, aby zachovala soucet T. Volba **Match total reco yield** nasledne normalizuje soucet cele reco predikce vcetne fakes na soucet dat. **Overall yield factor** je dalsi nasobek. **Fake normalization** je dodatecny faktor aplikovany az po prevazeni fake sablony.

V grafech je fialove aktualne prevazene fake spektrum. Stavovy radek uvadi aktualni fake fraction a nominalni podily fakes mimo truth rozsah / s nedostupnymi souradnicemi. Tyto podily se pocitaji z nominalnich vazenych MC poctu; nejsou to nove statisticke odhady. V 2D se podle aktivnich c,d vypisuji pT nebo spolecne truth kategorie.

Kazdy 1D/2D model se normalizuje podle vlastniho truth histogramu. Stejne parametry na ruznych modelech proto automaticky neznamenaji totozny globalni normalizacni faktor pro cely mnohorozmerny MC soubor. Pokud chcete prenest jednu 1D pT vahu do vsech response v Machine3.C, pouzijte konkretni ulozene 1D binove vahy.

## Sken a ukladani

**Scan prior grid...** funguje jako ve v2: konecna mrizka MC tilt, Power law nebo Both. Pouziva prave vybranou fake policy, stejny vypocet jako interaktivni GUI. V defaultu fake_factor i celkovy yield faktor zustavaji pevne. Minimum je minimum na zadane mrizce, nikoliv zaruka spojiteho globalniho minima.

Kriterium: `sum((pred-data)^2 / data_variance)` pres cely reco vektor, ve 2D pres vsechny bunky. Zaporne sWeighted obsahy se ponechaji. Biny s nulovou chybou se ze skore vylouci a jejich pocet se hlasi. Diagonalni skore nezahrnuje korelace sPlot ani nejistotu MC; neni to validovana p-hodnota.

**Save scan report** ulozi JSON a top-25 CSV vcetne fake policy a informaci o pokryti truth rozsahu. **Save current trial** ulozi JSON, NPZ (truth vahy, predikce i aktualni fake slozka) a PNG. **Load trial settings** obnovi i fake policy. Stare trial JSON bez fake_policy se interpretuji jako Fixed, aby se nezmenil jejich puvodni vyznam.

Pro prenos noveho vysledku do Machine3.C musi byt shodne nejen truth vahy, ale take zachazeni s Fake() a oblastmi mimo rozsah. Stary navod, ktery meni jen Fill/Miss a Fake nechava puvodni, reprodukuje pouze rezim Fixed. Novy nejlepsi prior je nutne znovu validovat v crossclosure.

## Prikazy

```bash
python3 prior_explorer_v3/prior_explorer_v3.py --prepare Unfolding/Machine.C --d0 1 5
python3 prior_explorer_v3/prior_explorer_v3.py PriorExplorerBundles/vas_novy_balicek.npz
python3 prior_explorer_v3/prior_explorer_v3.py --demo
python3 prior_explorer_v3/prior_explorer_v3.py --demo --scan --observable lambda1 --scan-family Both --fake-policy 'Truth weights / edge extension' --scan-output scan.json
python3 -m unittest discover -s prior_explorer_v3 -p 'test_prior_explorer*.py' -v
```

Pro cutBefore=false pridejte k priprave `--cut-after`, alternativni ROOT zadejte `--root /cesta/k/root`. CLI sken podporuje take `--centrality 0/1/2`, `--grid-json` s mapou `parametr: [min,max,krok]` a `--scan-output`. CLI vychazi z vychozich parametru; GUI z aktualnich posuvniku.

## Overeni

Proslo 30 numerickych testu, vcetne shody interaktivniho a davkoveho vypoctu pro vsechny fake policy a obe rodiny prioru, 1D/2D indexovani, konstantni response, signed data, kontinualniho limitu c,d -> 0 a explicitniho souctu fake kategorii.

C++ export byl zkompilovan a spusten v ROOT na syntetickem makru s minimalnim rozhranim RooUnfoldResponse. Zachyceni Fake(), vsech 24 modelu, underflow/overflow, chybejicich truth hodnot a prevod ROOT -> NPZ byly overeny. Nove fake predikce souhlasi s primym prevazenim jednotlivych syntetickych udalosti. Instrumentace byla take zkontrolovana na prilozenem Machine(7).C. Plny beh vaseho makra s vasimi vstupnimi MC/data soubory zde nebyl mozny.

Grafy byly vykresleny a vizualne zkontrolovany bez displeje. Zive Tk GUI v tomto prostredi otestovano nebylo.

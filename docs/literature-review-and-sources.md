# Изученная литература, протоколы и нормативные источники

Аннотированный каталог научных статей, лабораторных протоколов iGEM, патентов и российских отраслевых документов, использованных при проектировании конвейера DBTL и выборе целевой продуктовой платформы. Доступные в открытом доступе (Open Access) полнотекстовые PDF-версии ключевых статей сохранены в папках [`../papers/`](../papers/) и [`../Протоколы iGEM/`](../Протоколы%20iGEM/).

---

## 1. Локальные полнотекстовые статьи и протоколы в репозитории

### Биосинтез L-эрготионеина в дрожжах (ключевые референсы протокола)
1. **van der Hoek S. A., Darbani B., Zugaj K. E., Kell D. B., Borodina I. (2019).** *Engineering the yeast Saccharomyces cerevisiae for the production of L-(+)-ergothioneine.* Frontiers in Bioengineering and Biotechnology, 7:262. DOI: `10.3389/fbioe.2019.00262` (PMCID: `PMC6856651`, [Europe PMC](https://europepmc.org/article/PMC/PMC6856651)).
   * **Локальный полный текст:** [`papers/van_der_Hoek_2019_S_cerevisiae_ergothioneine_fulltext.txt`](../papers/van_der_Hoek_2019_S_cerevisiae_ergothioneine_fulltext.txt)
   * **Что взято в проект:** скрининг комбинаций грибных (`NcEgt1`, `CpEgt2`) и бактериальных (`MsEgtA–E`) генов в *S. cerevisiae* `CEN.PK113-7D`; доказательство того, что двухгенная связка *N. crassa Egt1 + C. purpurea Egt2* даёт 46–84 мг/л на базовом шаге и 598 мг/л в 1-л fed-batch реакторе; первичная проверка транспортера `AQR1` и водной экстракции кипячением.
2. **van der Hoek S. A., Rusnák M., Wang G., Jessop-Fabre M. M., Kell D. B., Borodina I. (2022).** *Engineering precursor supply for the high-level production of ergothioneine in Saccharomyces cerevisiae.* Metabolic Engineering, 70:129–142. DOI: `10.1016/j.ymben.2022.01.012` (PMCID: `PMC8967761`, [Europe PMC](https://europepmc.org/article/PMC/PMC8967761)).
   * **Локальный полный текст:** [`papers/van_der_Hoek_2022_S_cerevisiae_precursor_supply_2.39gL_fulltext.txt`](../papers/van_der_Hoek_2022_S_cerevisiae_precursor_supply_2.39gL_fulltext.txt)
   * **Что взято в проект:** точный состав минимальной среды Делфта без аминокислотных добавок, генетические мишени ветвей гистидина/SAM/цистеина (`HIS1`, `MET14`, `SAM2`, `SAH1`, делеция `SPE2`), достигнутый вручную потолок *de novo* на глюкозе **2,39 г/л** и зафиксированные авторами перекрёстные токсичности избытка серной ветви — прямая точка приложения для многомерной ML-модели DBTL.
3. **van der Hoek S. A., Rusnák M., Jacobsen I. H., Martínez J. L., Kell D. B., Borodina I. (2022).** *Engineering ergothioneine production in Yarrowia lipolytica.* FEBS Letters, 596(10):1356–1364. DOI: `10.1002/1873-3468.14239` (PMCID: `PMC9322481`, [Europe PMC](https://europepmc.org/article/PMC/PMC9322481)).
   * **Локальный полный текст:** [`papers/van_der_Hoek_2022_Y_lipolytica_ergothioneine_1.63gL_fulltext.txt`](../papers/van_der_Hoek_2022_Y_lipolytica_ergothioneine_1.63gL_fulltext.txt)
   * **Что взято в проект:** перенос тех же двух генов `NcEgt1 + CpEgt2` в промышленные олеагинные дрожжи *Y. lipolytica* (154 мг/л в малом объёме и **1,63 г/л** *de novo* на глюкозе в 1-л биореакторе при добавке фосфата).
4. **Команда iGEM HBUT-China (2022).** *CRISPR/Cas9 integration of NcEgt1 and CpEgt2 in S. cerevisiae CEN.PK113-5D.*
   * **Локальный PDF:** [`Протоколы iGEM/HBUT-China 2022 — эрготионеин.pdf`](../Протоколы%20iGEM/HBUT-China%202022%20—%20эрготионеин.pdf)
   * **Что взято в проект:** компактный трёхстраничный студенческий протокол ПЦР-сборки, LiAc/ПЭГ-трансформации `CEN.PK113-5D`, интеграции в межгенный локус `PDC1–STU2`, сброса плазмиды на 5-FOA и хроматографического подтверждения (+40,6% от второй копии кассеты).

### Методология цикла DBTL, машинное обучение и альтернативные кандидаты
5. **DBTL-предфинал (обзорная методологическая рукопись проекта, 2026).** *DBTL-цикл как основной подход проектирования и конструирования живых систем для повышения конкурентоспособности биотехнологических разработок в эпоху биоэкономики.*
   * **Локальный файл:** [`DBTL-предфинал.docx`](../DBTL-предфинал.docx)
   * **Содержание:** системная архитектура стадий Design, Build, Test и Learn; роль геномных моделей, комбинаторных библиотек, иерархических аналитических воронок, стандартов SBOL/FAIR, учёта отрицательных результатов, техноэкономической оценки (TEA) и контекста нацпроекта РФ «Технологическое обеспечение биоэкономики».
6. **Carbonell P. et al. (2018).** *An automated Design-Build-Test-Learn pipeline for enhanced microbial production of fine chemicals.* Communications Biology, 1:66. DOI: `10.1038/s42003-018-0076-9` (`PMC6123781`, [Europe PMC](https://europepmc.org/article/PMC/PMC6123781)).
   * **Локальный полный текст:** [`papers/Carbonell_2018_Automated_DBTL_pipeline_fulltext.txt`](../papers/Carbonell_2018_Automated_DBTL_pipeline_fulltext.txt)
7. **HamediRad M., Chao R., Weisberg S. et al. (2019).** *Towards a fully automated algorithm driven platform for biosystems design (BioAutomata).* Nature Communications, 10:5150. DOI: `10.1038/s41467-019-13189-z` (`PMC6853932`, [Europe PMC](https://europepmc.org/article/PMC/PMC6853932)).
   * **Локальный полный текст:** [`papers/HamediRad_2019_BioAutomata_lycopene_DBTL_fulltext.txt`](../papers/HamediRad_2019_BioAutomata_lycopene_DBTL_fulltext.txt)
   * **Значение:** классическая демонстрация байесовской оптимизации (Gaussian Process) биосинтеза ликопина в замкнутом цикле DBTL — подтверждение как работоспособности алгоритма на малых выборках, так и исчерпанности научной новизны самого ликопина.
8. **Radivojević T., Costello Z., Workman K., Garcia Martin H. (2020).** *A machine learning Automated Recommendation Tool (ART) for synthetic biology.* Nature Communications, 11:4879. DOI: `10.1038/s41467-020-18008-4` (`PMC7519088`, [Europe PMC](https://europepmc.org/article/PMC/PMC7519088)).
   * **Локальный полный текст:** [`papers/Radivojevic_2020_ART_machine_learning_DBTL_fulltext.txt`](../papers/Radivojevic_2020_ART_machine_learning_DBTL_fulltext.txt)
   * **Значение:** референсная архитектура вероятностного ансамбля ML для планирования следующих раундов метаболической инженерии по таблицам «вектор экспрессии → титр».
9. **Ma Y., Liu N., Greisen P. et al. (2022).** *Removal of lycopene substrate inhibition enables high carotenoid productivity in Yarrowia lipolytica.* Nature Communications, 13:572. DOI: `10.1038/s41467-022-28277-w` (`PMC8803881`, [Europe PMC](https://europepmc.org/article/PMC/PMC8803881)).
   * **Локальный полный текст:** [`papers/Ma_2022_Y_lipolytica_lycopene_17.6gL_fulltext.txt`](../papers/Ma_2022_Y_lipolytica_lycopene_17.6gL_fulltext.txt)
   * **Значение:** рекордные титры каротиноидов в *Y. lipolytica* (17,6 г/л ликопина и 39,5 г/л $\beta$-каротина), использованные при расчёте техноэкономического потолка каротиноидного направления.
10. **Chen S. et al. (2023).** *High-level production of pullulan (140.2 g/L) and its biosynthesis regulation in Aureobasidium pullulans BL06.* Frontiers in Bioengineering and Biotechnology, 11:1131875 (`PMC10467119`, [Europe PMC](https://europepmc.org/article/PMC/PMC10467119)).
    * **Локальный полный текст:** [`papers/Chen_2023_Aureobasidium_pullulans_140gL_fulltext.txt`](../papers/Chen_2023_Aureobasidium_pullulans_140gL_fulltext.txt)
    * **Значение:** экспериментальные данные по экстремальной вязкости культурального бульона пуллулана (до 6 270 мПа·с) и сопутствующему меланиногенезу.
11. **Команда iGEM WWU Münster (2022).** *MonChassis — Monoterpene production in Saccharomyces cerevisiae (71 стр. подробных дрожжевых протоколов).*
    * **Локальный файл:** [`Протоколы iGEM/WWU Münster 2022 — MonChassis.pdf`](../Протоколы%20iGEM/WWU%20Münster%202022%20—%20MonChassis.pdf)

---

## 2. Полный библиографический список проекта (43 источника обновлённого отчёта)

Полный пронумерованный список `[1]–[43]` (включая работы 2024–2026 гг.: *Tu & Stephanopoulos 2026* по ML-балансировке $\beta$-каротина, *De Biaggi & Ledesma-Amaro 2026* по модели `eciYali5-GEM`, *Zhou et al. 2025* по кроцетину в *Y. lipolytica*, *Moreno-Paz et al. 2024* по ML-DBTL *п*-кумаровой кислоты, *Ma et al. 2026 (BUCT)* — 3,05 г/л *de novo* в *S. cerevisiae*, *Hu et al. 2025* и *Chen et al. 2026* — 8,03 г/л в *Y. lipolytica*, *Yan et al. 2025* и *Zhang et al. 2026* — 7,2–12,16 г/л в *E. coli*, нотификацию *FDA GRN 001270* компании *GeneIII* от января 2026 г., *Ihsan & Cronan 2026* по обходному пути биотина в *PNAS*, методики *Carlsson 1974* и *Maurer / Seebeck 2020–2021*, патенты RU 2788417 C1, RU 2550202 C2 и RU 2769175 C1, Приказ Минсельхоза РФ № 392, Распоряжение Правительства РФ № 2777-р, аудит *RNC Pharma* 2026 г. и карточки ГРЛС Минздрава РФ по МНН «Адеметионин» на ноябрь 2026 г.) приведён в разделе 6 итогового отчёта [`ergothioneine-dbtl-report.pdf`](../ergothioneine-dbtl-report.pdf), а критический разбор — в [`critical-review-and-factcheck.md`](./critical-review-and-factcheck.md).

> *Уведомление по лицензиям открытых баз данных:* полнотекстовые открытые публикации получены через [Europe PMC](https://europepmc.org/) (условия использования отдельных статей определяются их лицензиями Creative Commons CC BY / CC BY-NC, указанными на титульных страницах PDF); рендеринг 2D химических структур выполнен с помощью библиотеки RDKit (BSD 3-Clause License).

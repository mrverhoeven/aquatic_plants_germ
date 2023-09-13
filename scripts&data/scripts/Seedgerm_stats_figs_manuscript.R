#'---
#' title: "Germinability and Viability of 4 macrophyte's seeds"
#' author: "Mike Verhoeven & Jonah Bacon"
#' output: 
#'    html_document:
#'       toc: true
#'       theme: default
#'       toc_depth: 2
#'       toc_float:
#'           collapsed: false
#'---

# header ------------------------------------------------------------------

#' # Preface
#' 
#' This code will pull in the datasets and analyze results of experiments to
#' evaluate the viability and germination of seeds from four species of aquatic
#' plants. The seeds were collected for use in a seed mix to be use for seeding
#' a set of experimental plots. These experiments aim to understand the number 
#' of sprouts we might expect to get from the seeds and how we might boost that 
#' number using pre-treatments.
#' 
#' Associated article: 
#' 
#'  Verhoeven, M.R., Bacon, J.A., and D.J. Larkin. Seed traits and germination treatments for dormancy break of four aquatic plant species. Target Journal: Aquatic Botany.
#' 
#' To-do:
#' Calculate "corrected germination" results
#' Table 1 update
#' clean up libraries
#' drop the time to germ analysis
#' 
#'  
#' ### Load libraries:
#'  
    # load libraries - needs check for necessity
    library(ggplot2)  
    library(data.table)
    library(MASS)
    library(GGally)
    library(survival)
    library(car) 
    library(emmeans)
    library(performance)
    library(KMsurv)
    library(rstanarm)  
    library(shinystan)
    library("bayesplot")
    library(projpred)
    library(sjPlot)
    library(sjlabelled)
    library(sjmisc)
    # devtools::install_github("kassambara/survminer", build_vignettes = FALSE)
    library(ggpubr)
    library(survminer)
    library(ggridges)
    library(gridExtra)

# load in data ------------------------------------------------------------
#' # Load and Review Data
#' ## Load data files
    
    #load germination dataset
    germ <- fread(file = "scripts&data/data/input/Seeds_germination.csv", nrows = 640) # uses nrows to make sure extraneous blank rows aren't loaded
    #load viability dataset
    viab <- fread(file = "scripts&data/data/input/Seeds_viability.csv")

# explore & clean up germination ------------------------------------------
#' ## Germination cleanup
#' 
#' This dataset includes germination data for 640 seeds that were part of
#' Jonah's lab experiments working to break dormancy of four species of seed.
#' The same covariates as in the viability dataset show up here, too (massed and
#' photographed to get those covariates)
#' 
#' ### Review and prep data      
    names(germ)
    str(germ)
    
    #change date formats:
    germ[ ,Date_chamber := as.IDate(Date_chamber, format = "%d-%b-%Y"), ]
    germ[ ,Date_germinated := as.IDate(Date_germinated, format = "%d-%b-%Y"), ]
    
    #make Treatment column values factors
    germ$Treatment <- as.factor(germ$Treatment)
    
    #subtract dates to get time to germ
    germ[ ,N_day_germ := Date_germinated - Date_chamber] #New column with days-to-germination
    
    #and add a column for yes or no germinated
    germ[ ,Germ_yn := as.logical(N_day_germ), ]
    
    #replace NAs with "FALSE"
    germ$Germ_yn <- replace(germ$Germ_yn, is.na(germ$Germ_yn), FALSE)
    
    # species pot_amp is actually pot_ill (based on Mike Verhoeven ID of grown-out adult of "B" seed's plant)
    germ[ , .N, Species]
    germ[Species == "Pot_amp", Species := "Pot_ill"]
    germ[ ,Species := as.factor(Species), ]
    germ[ ,N_day_germ := as.integer(N_day_germ) ,]
    
    # clean up treatment data - needs to be two sep columns for GA and SC
    germ[ , .N , .(Species, Treatment)]
    germ[Treatment == "SC+GA" | Treatment == "SC", scarify := TRUE , ]
    germ[ is.na(scarify) , scarify := FALSE , ]
    germ[Treatment == "SC+GA" | Treatment == "GA", gibberellic := TRUE , ]
    germ[ is.na(gibberellic) , gibberellic := FALSE , ]
    #check
    germ[ , .N , .(Treatment, gibberellic, scarify)]
    
    #l:w ratio
    germ[ ,lw_ratio := Major_axis/Minor_axis , ]
    
    #center and scale size pars within each species (why within and not across? Because we'll use these parameters for indiv species models--we do not analyze the whole gamut, instead we analyze each species independently)
    germ[Species == "Pot_ill"  , sc_lw := scale(lw_ratio)]
    germ[Species == "Pot_nat"  , sc_lw := scale(lw_ratio)]
    germ[Species == "Nup_var"  , sc_lw := scale(lw_ratio)]
    germ[Species == "Bra_sch"  , sc_lw := scale(lw_ratio)]
    
    germ[Species == "Pot_ill"  , sc_mass := scale(Mass_grams)]
    germ[Species == "Pot_nat"  , sc_mass := scale(Mass_grams)]
    germ[Species == "Nup_var"  , sc_mass := scale(Mass_grams)]
    germ[Species == "Bra_sch"  , sc_mass := scale(Mass_grams)]
    
    

# explore & clean up viability ------------------------------------------
#' ## Viability cleanup
#' 
#' This dataset includes one record for each of 168 seeds that were massed and
#' photographed before being cut open and stained with tetrazolium to test for 
#' viability. Jonah scored viability based on red-ness of the embryo. A set of 
#' negative controls were autoclaved before staining, and we determined that 
#' based on those neg controls any staining (score < 6 ) meant that the embryo
#' was viable (these seeds are not in the dataset). In the dataset we have two
#' "positive controls" which were harvested in August the day before the
#' viability assay. We determined that these were bad positive controls because
#' they had immature embryos.
#' 
#' ### Review and prep data       

    names(viab)
    str(viab)
    
    # add a column for viable yes or no
    viab[ ,Viab_yn := as.logical(Embryo_stained < 6), ]
    
    # species pot_amp is actually pot_ill (based on Mike V ID of grown-out "B" seed's plant)
    viab[ , .N, Species]
    viab[Species == "Pot_amp", Species := "Pot_ill"]
    viab[Species == "pot_amp_cont", Species := "Pot_ill_cont"]
    viab[ ,Species := as.factor(Species), ]
    
    #l:w ratio
    viab[ ,lw_ratio := Major_axis/Minor_axis , ]
    
    #center and scale size pars (again done by species to fit with species by species modeling):
    viab[Species == "Pot_ill"  , sc_lw := scale(lw_ratio)]
    viab[Species == "Pot_nat"  , sc_lw := scale(lw_ratio)]
    viab[Species == "Nup_var"  , sc_lw := scale(lw_ratio)]
    viab[Species == "Bra_sch"  , sc_lw := scale(lw_ratio)]
    
    viab[Species == "Pot_ill"  , sc_mass := scale(Mass_grams)]
    viab[Species == "Pot_nat"  , sc_mass := scale(Mass_grams)]
    viab[Species == "Nup_var"  , sc_mass := scale(Mass_grams)]
    viab[Species == "Bra_sch"  , sc_mass := scale(Mass_grams)]
    
    #drop the data for the visual pos controls (these were used to set the scale for our TZ stain, but are not to be used for analysis):
    viab <- viab[!Species %in% c('Pot_ill_cont', 'Nup_var_cont')]
    viab[ , Species := droplevels(Species)]
    

# viability ---------------------------------------------------------------

#' # Analysis of Viability
#' ## Visualize
#' 
    #prep fields for figure:
    viab[is.na(Viab_yn) , viab_words := "No embryo" ]
    viab[Viab_yn == "TRUE" , viab_words := "Viable" ]
    viab[Viab_yn == "FALSE" , viab_words := "Non-viable"  ]
        viab[, summary(viab_words)]
    
    #set ggplot theme 
    theme_set(theme_bw(base_size = 14))
    
    #write plot
    viabilityplot <- ggplot(viab, aes(x = Species, fill = viab_words))+
        geom_bar(position = "stack") +
        coord_flip()+
        scale_fill_manual(values=c("#D55E00", "#E69F00", "#009E73")) +
        #geom_text() +
        geom_text(aes(label = after_stat(count)), stat = "count", size = 5, position = position_stack(vjust = 0.5), colour = "black")+
        labs(fill = NULL)  +
        scale_y_continuous(expand = c(0,0)) +
        ylab("Number of Seeds") +
        xlab("")+
        theme_bw() +
        scale_x_discrete(position = "top", limits=c("Bra_sch", "Nup_var", "Pot_ill", "Pot_nat" ),
                         labels = c("     Brasenia schreberi", "       Nuphar variegata", "   Potamogeton illinoiensis", "    Potamogeton natans"))+
        theme(
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            text = element_text(size = 14, color = "black"),
            axis.text.y = element_text(size = 12, face = "italic", hjust = 1),
            legend.position = c(0.85,0.4),
            legend.background = element_blank())
    
    #view plot
    viabilityplot
    
# bayesian viability ------------------------------------------------------

#' ## Model viability results
#' Bayesian estimation of posteriors for viability (resp ~ 0 + pred will suppress intercepts, resp ~ pred - 1 gives means coded model matrix) using binomial error distributions. Results are log-odds, and use back-transformation (plogis) to get back to the % viable scale:
    
    #viability in bayesian framework:
    
    #reorder species list
    viab[, Species := factor(Species, levels = c('Pot_nat', 'Pot_ill', 'Nup_var', 'Bra_sch'))]
    
    #write model
    bayes_species_viab <- stan_glm(Viab_yn ~ Species - 1  , data = viab[Species %in% c('Pot_ill','Bra_sch','Nup_var','Pot_nat'), ,], family = binomial(),  set.seed(4913))
    summary(bayes_species_viab)
    
    #results in fig
    viabmodelfig <- plot_model(bayes_species_viab, 
                               bpe = "median", 
                               prob.inner = 0.5, 
                               prob.outer = 0.95, 
                               transform = 'plogis', 
                               show.values = T,
                               title = "", axis.title = c("Proportion viable", NULL), 
                               axis.lim = c(0,1), 
                               axis.labels = c(" ", " ", " "," "),
                               colors = "bw",
                               text = element_text(size = 14, color = "black"))
    #view figure
    viabmodelfig 
    

    
    #print model results to table:                      
    ci95 <- plogis(posterior_interval(bayes_species_viab, prob = 0.95))
    bayes_species_viab_ci_tab <- data.table(round(ci95, 2))
    bayes_species_viab_ci_tab[ , species := c("Potamogeton natans","Potamogeton illinoensis", "Nuphar variegata", "Brasenia schreberi"  ) , ]
    bayes_species_viab_ci_tab[ , mean := round(plogis(bayes_species_viab$coefficients),2) , ]
    
    bayes_species_viab_ci_tab[, maxgerm := round(c(11/40, 14/40, 1/40, 0/40),2)] #GA and Scarif rates
    bayes_species_viab_ci_tab[, germ_of_viab := round(maxgerm/mean, 2)]
    
    #view simple results output
    bayes_species_viab_ci_tab
    
    # Write viability results to table
    #fwrite(bayes_species_viab_ci_tab, file = "viability_table.csv")
    
    #create figure 1:

# Viability Viz: ---------------------------------------------------------------

    #model diagnostics:
#' ## Bayesian Viability Model Diagnostics
#' We follow https://mc-stan.org/users/documentation/case-studies/rstan_workflow.html for diagnostics:

# launch_shinystan(bayes_species_viab)

#' ## Interpret viability data by species:
#' 
#' Plot designed following: https://cran.r-project.org/web/packages/sjPlot/vignettes/plot_model_estimates.html
#' We can see that all species had non-zero viability. We can use 
#' these median value for each species to discount our germ in the 
#' germination studies. We will call this a "realized" germination.
#' Basically, because only 53% of the B screberi seeds were viable, we only 
#' really tested germination in 53% of the 40 seeds because 47% were
#' nonviable anyway and could thus NEVER germinate. 
#' 
#' ## Viability Viz                
    #two plot components:
    viabmodelfig 
    viabilityplot
    
    Viabviz <- ggarrange( viabilityplot, ggplot() + theme_void(), viabmodelfig,
                          widths = c(2,-0.06, 1),
                          labels = c("A", "", "B"),
                          label.x = c(0.01,0,0.05),
                          label.y = c(0.99,0,0.99),
                          font.label = list(size = 18, color = "black", face = "bold", family = NULL),
                          nrow = 1)
    Viabviz
    
    # write to file:
    # # Customizing the output
    # tiff(filename = "Figure_1.tiff",         # File name
    #      width = 11, height = 4, units = "in", # Width and height in inches
    #      # bg = "white",          # Background color
    #      # colormodel = "cmyk",   # Color model (cmyk is required for most publications)
    #      # paper = "A4"
    #      # pointsize = 12,
    #      # compression = c("none", "rle", "lzw", "jpeg", "zip", "lzw+p", "zip+p"),
    #      res = 300, family = "", restoreConsole = TRUE,
    #      type =  "cairo"
    # )
    # 
    # # Creating a plot
    # Figure1
    # 
    # # Closing the graphical device
    # dev.off()




    


# analysis of germ --------------------------------------------------------

#' # Analysis of Germination
#' ## Uncorrected Germ Data
#' NOt corrected for viability 
    
    # Code for stacked barplot
    germ[, Species := factor(Species, levels = c('Pot_nat', 'Pot_ill',  'Nup_var', 'Bra_sch'))]
    species.names <- c("P. natans","P. illinoensis","N. variegata","B. schreberi" )
    names(species.names) <- c("Pot_nat","Pot_ill","Nup_var","Bra_sch"  )
    
    #write out plot
    uncorr_germ_viz <- ggplot(germ[ , .N, .(Species, Germ_yn, Treatment)] , aes(x = Treatment, y = N , fill = Germ_yn, label = N)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(values=c("#D55E00",  "#009E73"), labels = c("Ungerminated", "Germinated")) +
        facet_grid(vars(group = Species), labeller = labeller(group = species.names)) +
        geom_text(size = 5, position = position_stack(vjust = 0.85)) +
        labs(fill = "Germinated?") +
        scale_x_discrete( breaks=c("CT","GA","SC","SC+GA"), labels=c("Control", "GA", "Scarified", "Scarified+ \n GA")) +
        ylab("Number of Seeds") +
        theme(
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            strip.text.y = element_text(size = 12, face = "bold.italic"),
            panel.spacing = unit(1, "lines"),
            text = element_text(size = 12, face = "bold", color = "black"),
            legend.title = element_blank()
        )
    
    #view plot
    uncorr_germ_viz


    # # Write fig to file
    # # Customizing the output
    # tiff(filename = "Figure_2.tiff",         # File name
    #      width = 5.2, height = 6.5, units = "in", # Width and height in inches
    #      # bg = "white",          # Background color
    #      # colormodel = "cmyk",   # Color model (cmyk is required for most publications)
    #      # paper = "A4"
    #      # pointsize = 12,
    #      # compression = c("none", "rle", "lzw", "jpeg", "zip", "lzw+p", "zip+p"),
    #      res = 300, family = "", restoreConsole = TRUE,
    #      type =  "cairo"
    # )
    # 
    # # Creating a plot
    # Figure2
    # 
    # # Closing the graphical device
    # dev.off()

#' ## How to analyze?
#' 
#' We can see two important things from this figure that we want to put stat
#' values on:
#' 
#'  1.  Scarification and gibberellic acid look like they had a Synergistic
#'  interaction (e.g. in P. natans the effect of scarifying and GA treatments 14
#'  was more than the effect of adding the result of scarifiction 0 and GA 5)
#'  2. The effect of treatments & traits varied across species, with B Schreberi and N
#'  variegata looking similar, and both Potamogetons looking similar.
#'  3. We need to account for the viability in the germination results-- in practice, 
#'  that means dividing the germ by the viability number. We could try a hierarchical model 
#'  to assign the non-germ variance across both levels, but we really do not have sample sizes 
#'  for that sorta operation. Instead, let's estimate everything as previously done, then 
#'  do a post-hoc correction for viability.  
#'  
#'  Thus, for each species, we can evaluate these treatment hypotheses using a 
#'  model for germination. IN our reporting we will use viability to generate
#'  corrected germination.   
#'  
#'  Results for binomial GLMs in the frequentist OLS pumped out odd results. 
#'  Digging into the seed germ analysis lit suggests that the high occurrence of
#'  0% germ in our data (termed "complete separation") means that there's no
#'  variance in those categories for a
#'  model to fit to... And the minimum sample size of n*p > 5 for a binomial
#'  estimator tells us that for a 1% germination rate we should have used a 
#'  sample size of 500. Because we saw at best 1/1185 Pot nat sprouting
#'  (sprouting experiment with 0 sprouts), we should have run a sample size of
#'  at least 5925 (5/(1/1185)) for the Pot nat controls (no GA or Scarif)....
#'  
#'  So, with all of that in consideration, the hyp testing and the model 
#'  selection is prob inappropriate with the data that we collected being
#'  evaluated with a binomial glm in a frequentist framework.See: 
#'  
#'  Gianinetti, A. (2020). Basic features of the analysis of germination
#'   data with generalized linear mixed models. Data, 5(1). 
#'   https://doi.org/10.3390/data5010006
#'   
#'  McNair, J. N., Sunkara, A., & Frobish, D. (2012). How to analyse seed 
#'   germination data using statistical time-to-event analysis: Non-parametric
#'   and semi-parametric methods. Seed Science Research, Vol. 22, pp. 77–95. 
#'   https://doi.org/10.1017/S0960258511000547
#'  
#'  The issue at hand is that we have categories that produced true zeroes.
#'  While this wasn't entirely unexpected, it does result in a model that has no
#'  real variance estimate for the zero categories (SE estimates of Infinity for
#'  those categories).
#'  
#'  One way for us to overcome this is to avoid the issue of calculating the
#'  variance using frequentist approaches and switch to estimating the models 
#'  via a bayesian implementation of a glm.  
#'  
#'  More info: 
#'  https://stats.stackexchange.com/questions/381915/glm-for-count-data-with-all-zeroes-in-one-category
#'  
#'  rstanrm vignette:
#'  http://mc-stan.org/rstanarm/articles/rstanarm.html
#'  
#'  rstanrm user guide:
#'  http://dx.doi.org/10.20982/tqmp.14.2.p099
#'  
#' 


# hypothesis testing by species ---------------------------------------------
#' ## Potamogeton natans 
    #Pot nat:
    #model
    bayes_pnat <- stan_glm(data = germ[Species == "Pot_nat"],  Germ_yn ~ gibberellic * scarify + sc_mass + sc_lw, family = binomial(link = 'logit'), set.seed(4913))
    #output
    summary(bayes_pnat)
    #visualize results on data scale
    plot_model(bayes_pnat, prob = 0.5, prob_outer = 0.95, transform = 'plogis', axis.lim = c(0,1), title = "Potamogeton natans")
    #visualize on model scale (log-odds)
    plot_model(bayes_pnat, prob = 0.5, prob_outer = 0.95, title = "Potamogeton natans", colors = "black", transform = NULL,
               axis.labels = c("GA3 + Scarification", "Length:Width", "Mass", "Scarify", "GA3"))
    
    #print model results to table:                      
    ci.table <- cbind(posterior_interval(bayes_pnat, prob = 0.50),posterior_interval(bayes_pnat, prob = 0.95))
    bayes_germ_ci_tab <- data.table(species = "Potamogeton natans",
                                    parameter = rownames(ci.table),
                                    est = bayes_pnat$coefficients,
                                    ci.table
                                    )
    bayes_germ_ci_tab_master <- bayes_germ_ci_tab

#' ## Potamogeton illinoensis
    #Pot ill
    bayes_pill <- stan_glm(data = germ[Species == "Pot_ill"], Germ_yn ~ gibberellic * scarify + sc_mass + sc_lw  ,  family = binomial(link = 'logit'), set.seed(4913))
    summary(bayes_pill)
    #visualize model results
    plot_model(bayes_pill, title = "Potamogeton illinoensis", colors = "black", transform = NULL)
    plot_model(bayes_pill, title = "Potamogeton illinoensis", colors = "black", transform = "exp", ci.style = "bar")
    plot_model(bayes_pill, prob = 0.5, prob_outer = 0.95, transform = 'plogis', axis.lim = c(0,1), title = "Potamogeton illinoensis",vline = 0.5, vline.color = "black")
    
    #print model results to table:
    ci.table <- cbind(posterior_interval(bayes_pill, prob = 0.50),posterior_interval(bayes_pill, prob = 0.95))
    bayes_germ_ci_tab <- data.table(species = "Potamogeton illinoensis",
                                    parameter = rownames(ci.table),
                                    est = bayes_pill$coefficients,
                                    ci.table
    )
    
    bayes_germ_ci_tab_master <- rbind(bayes_germ_ci_tab_master, bayes_germ_ci_tab)

#' ## Nuphar variegata 
#' we have decided that while we *can* estimate the effects for N var, our lack
#' of overall germ should be used to indicate that it's not a real good idea to
#' use these data to estimate trt or trait effects..
        
    #nup var
    #we have decided that while we can estimate the efects for N var, our lack of overall germ should be 
    # bayes_nvar <- stan_glm(Germ_yn ~ gibberellic * scarify + sc_mass + sc_lw, data = germ[Species == "Nup_var"], family = binomial(link = 'logit'), set.seed(4913))
    # summary(bayes_nvar) # no sig trt effects
    # plot(bayes_nvar)
    # plot_model(bayes_nvar, type = 'est', prob = 0.5, prob_outer = 0.5,
    #            title = "Nuphar variegata",
    #            colors = "black", transform = NULL)
    # plot_model(bayes_nvar, prob = 0.5, prob_outer = 0.95, transform = 'plogis', axis.lim = c(0,1), title = "Nuphar variegata")
    # 
    # 
    # #print model results to table:
    # ci.table <- cbind(posterior_interval(bayes_nvar, prob = 0.50),posterior_interval(bayes_nvar, prob = 0.95))
    # bayes_germ_ci_tab <- data.table(species = "Nuphar variegata",
    #                                 parameter = rownames(ci95),
    #                                 est = bayes_nvar$coefficients,
    #                                 ci.table
    # )
    # 
    # bayes_germ_ci_tab_master <- rbind(bayes_germ_ci_tab_master, bayes_germ_ci_tab)


    bayes_germ_ci_tab_master[, exp(est),]
    
#' ## Brasenia schreberi 
#' Because Brasenia sample includes no germ in any treatment or set of traits, 
#' this model will not run. There's nothing to estimate on.  
    
    # bra schr
    # Brasenia sample includes no germ--this model will not run
    # bayes_bsch <- stan_glm(Germ_yn ~ gibberellic * scarify + sc_mass + sc_lw , data = germ[Species == "Bra_sch"], family = binomial(link = 'logit'), set.seed(4913)) # all zeros = not gonna work
    
    # Write CI table to file
    # write.csv(bayes_germ_ci_tab_master, file = "bayes_germ_ci_table.csv")



# create Figure 2 ---------------------------------------------------------

#' ## Figure 2    
    
bayes_germ_ci_tab_master[ , parameter := factor(parameter, levels = c(rev(unique(bayes_germ_ci_tab_master$parameter)))) ,]

    bayes_germ_ci_tab_master[ , species := factor(species, levels = c(unique(bayes_germ_ci_tab_master$species))) ,]
    bayes_germ_ci_tab_master[ parameter == "gibberellicTRUE:scarifyTRUE", letter := c(" A", "B") ]
    
    Figure2 <- ggplot(bayes_germ_ci_tab_master[parameter != "(Intercept)"], aes( parameter, est ))+
        geom_crossbar(aes(ymax = bayes_germ_ci_tab_master[parameter != "(Intercept)"]$'75%',
                          ymin = bayes_germ_ci_tab_master[parameter != "(Intercept)"]$'25%'),
                      fill = "gray", 
                      fatten = 1.5, 
                      width = .5)+
        geom_linerange(aes(ymax = bayes_germ_ci_tab_master[parameter != "(Intercept)"]$'97.5%',
                           ymin = bayes_germ_ci_tab_master[parameter != "(Intercept)"]$'2.5%'), orientation = "x")+
        facet_wrap(~species)+
        geom_vline(xintercept = 0)+
        geom_hline(yintercept = 0)+
        ylab("Log-Odds")+
        xlab(NULL)+
        scale_x_discrete(labels=c("GA3:Scarification", "Length/Width Ratio", "Mass", "Scarify", "GA3", "Intercept"))+
        theme(axis.text.x =  element_text())+
        theme(
            strip.text.x = element_text(
                size = 12, face = "italic"
            ))+
        theme(
            strip.background = element_rect(
                color= "black", fill= "white", size=.5, linetype="solid"
            )
        )+
        ylim(c(-13,10))+
        # geom_text(aes(label = letter), vjust = -0.1, hjust = 3, size = 15)+
        coord_flip()
    Figure2
    
    # # Customizing the output
    tiff(filename = "Figure_2.tiff",         # File name
         width = 7, height = 3.25, units = "in", # Width and height in inches
         # bg = "white",          # Background color
         # colormodel = "cmyk",   # Color model (cmyk is required for most publications)
         # paper = "A4"
         # pointsize = 12,
         compression = c("lzw"),
         res = 300, family = "", restoreConsole = TRUE,
         type =  "cairo"
    )

    # Creating a plot
    Figure2

    # Closing the graphical device
    dev.off()
    
    
# bayesian model diagnostics ----------------------------------------------


#' 
#' ## Bayesian Germinability Model Diagnostics
#' We follow 
#' https://mc-stan.org/users/documentation/case-studies/rstan_workflow.html for
#' diagnostics:
    
    # launch_shinystan(bayes_pnat)
    # launch_shinystan(bayes_pill)
    # launch_shinystan(bayes_nvar)


# trait distributions -----------------------------------------------------


        
    ggplot(germ, aes(x = Mass_grams, y = Species)) + geom_density_ridges()
    
    ggplot(germ[Area < 4000000], aes(x = Area, y = Species)) + geom_density_ridges() 
    
    ggplot(germ, aes(x = lw_ratio, y = Species)) + geom_density_ridges()
    
    ggplot(germ, aes(x = Red_mean, y = Species)) + geom_density_ridges() 
    
    
    names(germ)[c(18:22,25,28,36)]
        
       plot_dat <-  melt(germ[Area< 4000000 & Green_mean < 120 , ,], id.vars = "Species", measure.vars = names(germ)[c(18:22,25,28,36)], variable.name = "trait", value.name = "value") 
       
       plot_dat <- plot_dat[trait %in% c("Mass_grams", "lw_ratio", "Major_axis", "Minor_axis")]
       
       plot_dat[, trait := factor(trait, levels = c("Mass_grams", "lw_ratio", "Major_axis", "Minor_axis")) ]
       
       plot_dat[trait == "Mass_grams" , trait := "Mass (g)"]
       plot_dat[trait == "lw_ratio" , trait := "Length/Width Ratio"]
       plot_dat[trait %in% c("Major_axis", "Minor_axis"), value := value/1000]
       plot_dat[trait == "Major_axis" , trait := "Length (mm)"]
       plot_dat[trait == "Minor_axis" , trait := "Width (mm)"]
       

    Figure1 <- ggplot(plot_dat[], aes(value, Species))+
        geom_density_ridges(aes( fill =  Species))+
        facet_wrap(~trait, scales = "free" , ncol = 2)+
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              axis.ticks.y=element_blank(),
              axis.text.y = element_blank(),
              legend.text=element_text(face = "italic"))+
        xlab("")+
        ylab("")+
        scale_fill_discrete(position = "top", limits=c("Bra_sch", "Nup_var", "Pot_ill", "Pot_nat" ),
                         labels = c("Brasenia schreberi", "Nuphar variegata", "Potamogeton illinoiensis", "Potamogeton natans"))
    
    Figure1
    
    # # Customizing the output
    tiff(filename = "Figure_1.tiff",         # File name
         width = 6, height = 5, units = "in", # Width and height in inches
         # bg = "white",          # Background color
         # colormodel = "cmyk",   # Color model (cmyk is required for most publications)
         # paper = "A4"
         # pointsize = 12,
         compression = c("lzw"),
         res = 300, family = "", restoreConsole = TRUE,
         type =  "cairo"
    )
    
    # Creating a plot
    Figure1
    
    # Closing the graphical device
    dev.off()
    
    
    
    
        


# Table 1: relative germinability -----------------------------------------

    #viability in bayesian framework:
    
    #reorder species list
    germ[ , .N , Species]
        germ[, Species := factor(Species, levels = c('Pot_nat', 'Pot_ill', 'Nup_var', 'Bra_sch'))]
    
    #write model
    bayes_species_germ <- stan_glm(Germ_yn ~ Species - 1  , data = germ[Species %in% c('Pot_ill','Pot_nat') & Treatment == "SC+GA", ,], family = binomial(),  set.seed(4913))
    summary(bayes_species_germ)
    
    #results in fig
   plot_model(bayes_species_germ, 
                               bpe = "median", 
                               prob.inner = 0.5, 
                               prob.outer = 0.95, 
                               transform = 'plogis', 
                               show.values = T,
                               title = "", axis.title = c("Proportion viable", NULL), 
                               axis.lim = c(0,1), 
                               # axis.labels = c(" ", " ", " "," "),
                               colors = "bw",
                               text = element_text(size = 14, color = "black"))

    
    
    
    #print model results to table:                      
    ci95 <- plogis(posterior_interval(bayes_species_germ, prob = 0.95))
    bayes_2species_germ_ci_tab <- data.table(round(ci95, 2))
    bayes_2species_germ_ci_tab[ , species := c("Potamogeton natans","Potamogeton illinoensis") , ]
    bayes_2species_germ_ci_tab[ , mean := round(plogis(bayes_species_germ$coefficients),2) , ]
    
    bayes_2species_germ_ci_tab[, maxgerm := round(c(11/40, 14/40),2)] #GA and Scarif rates
    bayes_2species_germ_ci_tab[, germ_of_viab := round(maxgerm/mean, 2)]
    
    #view simple results output
    bayes_species_viab_ci_tab    
    
    
# footer ------------------------------------------------------------------
#' ## Document footer 
#' 
#' Session Information:
#+ sessionInfo
sessionInfo()

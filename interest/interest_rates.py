import numpy as np


def ISA_calculations(invest_month,years):
    # returns total_value
    # outputs total value, net_contributions, surplus
    interest_rate = 0.07
    government_contribution = 0.25

    invest_LISA = invest_month / 2
    invest_ISA = invest_month / 2


    principle_ISA = 3000
    principle_LISA = 0
    total_value = principle_ISA + principle_LISA

    net_contributions = 0

    months = int(years * 12)

    for i in range(months):
        principle_ISA = (principle_ISA + invest_ISA) * (1 + interest_rate/12)
        principle_LISA = (principle_LISA + invest_LISA * (1+government_contribution)) * (1 + interest_rate/12)
        total_value = principle_ISA + principle_LISA
        net_contributions += invest_month 


    print(total_value)
    print(net_contributions)
    print('surplus : ',total_value-net_contributions)
    
    return total_value

def income_calculations():
    GTA_rate = 17.75
    GTA_hours_semester = 36
    GTA_yearly = GTA_rate * GTA_hours_semester
    GTA_monthly = GTA_yearly / 12

    STIPEND_monthly = 1603.08

    EXTERNAL_monthly = 500

    TOTAL_month = GTA_monthly + STIPEND_monthly + EXTERNAL_monthly
    print('income per month : ',TOTAL_month)
    print('income per year : ',12*TOTAL_month)
    return TOTAL_month

def outcome_calculations(INCOME):
    def hobbies():
        def ice_skating():
            PLANET_ICE_FEES_WEEK = 10 * 3
            COACHING_FEES_WEEK = 25 
            SKATES_YEAR = 200

            total_month = SKATES_YEAR / 12 + PLANET_ICE_FEES_WEEK * 4.333333 + COACHING_FEES_WEEK * 4.333333
            return total_month
        def archery():
            bruntwood_year = 100
            uni_year = 100

            total_month = (bruntwood_year + uni_year) / 12  
            return total_month
        
        total_hobbies = ice_skating() + archery()

        return total_hobbies
    
    def utilities_rent():
        rent_month = 960
        utilities = 100 + 50

        total_month = rent_month + utilities
        return total_month

    def eating():
        week = 50
        total_month = 50*4.3333333
        return total_month


    total_outcome_month = hobbies() + utilities_rent() + eating() + 500
        
    print(total_outcome_month)


INCOME = income_calculations()
ISA_value = ISA_calculations(500,3)

outcome_calculations(INCOME)
import numpy as np


def ISA_calculations(initial_ISA,initial_LISA,invest_month,years):
    # returns total_value
    # outputs total value, net_contributions, surplus
    interest_rate = 0.1
    government_contribution = 0.25

    invest_LISA = invest_month / 2 
    invest_ISA = invest_month / 2


    principle_ISA = initial_ISA
    principle_LISA = initial_LISA
    total_value = principle_ISA + principle_LISA

    net_contributions = 0

    months = int(years * 12)

    for i in range(months):
        principle_ISA = (principle_ISA + invest_ISA) * (1 + interest_rate/12)
        principle_LISA = (principle_LISA + invest_LISA * (1+government_contribution)) * (1 + interest_rate/12)
        total_value = principle_ISA + principle_LISA
        net_contributions += invest_month 


    # print(total_value)
    # print(net_contributions)
    # print('surplus : ',total_value-net_contributions)
    
    return total_value

def income_calculations():
    GTA_rate = 17.75
    GTA_hours_semester = 36
    GTA_yearly = GTA_rate * GTA_hours_semester
    GTA_monthly = GTA_yearly / 12

    STIPEND_monthly = 1603.08

    EXTERNAL_monthly = 500

    ADDITIONAL_AI_weekly = 100
    ADDITIONAL_AI_monthly = ADDITIONAL_AI_weekly * 4.3333333

    TOTAL_month = GTA_monthly + STIPEND_monthly + EXTERNAL_monthly + ADDITIONAL_AI_monthly
    # print('income per month : ',TOTAL_month)
    # print('income per year : ',12*TOTAL_month)
    return TOTAL_month

def outcome_calculations(INCOME,invest):
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
        virgin = 33.74
        gas_electric = 65.70
        water = 30 # tentative
        utilities = virgin + gas_electric + water
        total_month = rent_month + utilities
        return total_month

    def eating():
        fried_rice = (1.70 + 0.65 + (2.39/3*2) + 1.79/5)/2
        pizza = 1.99 + 0.74/2

        meal_deal = 3.6
        
        sandwiches = 0.99/2 + 3/7


        week = fried_rice * 6 + pizza + sandwiches * 7


        total_month = week*4.3333333
        return total_month

    total_outcome_month = hobbies() + utilities_rent() + eating() + invest
        
    
    return total_outcome_month

INCOME = income_calculations()
invest = 1100

ISA_value = ISA_calculations(3000,0,invest,40) # 3 is length of PhD 
OUTCOME = outcome_calculations(INCOME,invest)

print('net income : ', INCOME)
print('net outgoings : ',OUTCOME)
print('net difference : ',INCOME - OUTCOME)
print('Estimated ISA value at end of PhD : ',ISA_value)
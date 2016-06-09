################################################################################
# Nicolas Tancogne-Dejean, 2012-2013
################################################################################

RM := rm -rf

# All Target
all: doc

# Tool invocations
doc:
        @echo "Generating Documentation"
        doxygen doxygen.in > doxygen.log 

# Other Targets
clean:
        $(RM) html latex rtf *.log
        @echo ' '

.PHONY: clean
~                    

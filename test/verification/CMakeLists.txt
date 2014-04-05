#
# Add the tests
#
add_test (Dakota_2x1x1_np1 python ./Dakota_2x1x1/np1/Dakota_2x1x1.py)
add_test (ElasticPlasticCorrespondenceFullyPrescribedTension_np1 python ./ElasticPlasticCorrespondenceFullyPrescribedTension/np1/ElasticPlasticCorrespondenceFullyPrescribedTension.py)

add_custom_target( vtest
   COMMAND ctest
)

# Assembly note

All 18 mechanical integrations completed successfully. The first evidence assembly used a nonexistent `Mixture.energy` alias. It failed after writing the mechanical results; the call was corrected to the existing analytic `Mixture.du` API and the complete evidence assembly was rerun. No mechanics, coefficients or accepted data were changed.

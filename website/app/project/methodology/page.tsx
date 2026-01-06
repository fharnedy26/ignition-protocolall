import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";

export default function MethodologyPage() {
  return (
    <div className="container mx-auto px-4 py-12">
      <div className="mb-12">
        <h1 className="mb-4 text-5xl font-bold">
          <span className="bg-gradient-to-r from-primary to-accent bg-clip-text text-transparent">
            Methodology
          </span>
        </h1>
        <p className="text-xl text-muted-foreground">
          High-level overview of the Ignition Protocol approach to fuel discovery and optimization.
        </p>
      </div>

      <div className="space-y-8">
        {/* System Architecture */}
        <Card>
          <CardHeader>
            <CardTitle className="text-2xl">System Architecture</CardTitle>
          </CardHeader>
          <CardContent className="space-y-4">
            <p className="text-muted-foreground">
              Ignition Protocol is a unified pipeline with two coupled engines:
            </p>
            <div className="grid gap-4 md:grid-cols-2">
              <div className="rounded-lg border border-border bg-card p-4">
                <h3 className="mb-2 font-semibold text-primary">Molecular Generation</h3>
                <ul className="space-y-1 text-sm text-muted-foreground">
                  <li>• Fragment-based evolution</li>
                  <li>• RDKit validation</li>
                  <li>• Multi-objective fitness</li>
                  <li>• Diversity preservation</li>
                </ul>
              </div>
              <div className="rounded-lg border border-border bg-card p-4">
                <h3 className="mb-2 font-semibold text-primary">Blend Optimization</h3>
                <ul className="space-y-1 text-sm text-muted-foreground">
                  <li>• Deterministic algorithm</li>
                  <li>• Fuel-type constraints</li>
                  <li>• Lab-ready outputs</li>
                  <li>• Reproducible results</li>
                </ul>
              </div>
            </div>
          </CardContent>
        </Card>

        {/* Molecular Generation */}
        <Card>
          <CardHeader>
            <CardTitle className="text-2xl">Molecular Generation Module</CardTitle>
          </CardHeader>
          <CardContent className="space-y-4">
            <div>
              <h3 className="mb-2 font-semibold text-primary">Fragment-Based Evolution</h3>
              <p className="text-muted-foreground">
                The system uses fragment libraries (GDB-11 database) to evolve novel fuel molecules. Fragments 
                are combined using genetic algorithms with crossover, mutation, and selection operators.
              </p>
            </div>
            <div>
              <h3 className="mb-2 font-semibold text-primary">Fitness Function</h3>
              <p className="text-muted-foreground">
                Candidates are evaluated through a multi-objective fitness function that balances:
              </p>
              <ul className="mt-2 space-y-1 text-muted-foreground">
                <li>• Molecular weight</li>
                <li>• Oxygen-to-carbon ratio</li>
                <li>• Volatility (logP)</li>
                <li>• Combustion potential</li>
                <li>• H-bonding capacity</li>
              </ul>
            </div>
            <div>
              <h3 className="mb-2 font-semibold text-primary">Validation & Filtering</h3>
              <p className="text-muted-foreground">
                All molecules are validated using RDKit to ensure 100% chemical validity. The system filters for 
                fuel-relevance (C/H/O content) and preserves diversity through Tanimoto similarity metrics.
              </p>
            </div>
          </CardContent>
        </Card>

        {/* Blend Optimization */}
        <Card>
          <CardHeader>
            <CardTitle className="text-2xl">Blend Optimization Engine</CardTitle>
          </CardHeader>
          <CardContent className="space-y-4">
            <div>
              <h3 className="mb-2 font-semibold text-primary">Greedy-Then-Refine Algorithm</h3>
              <p className="text-muted-foreground">
                The optimizer uses a deterministic greedy selection followed by refinement. The algorithm:
              </p>
              <ul className="mt-2 space-y-1 text-muted-foreground">
                <li>• Selects components using multi-objective scoring (RON, cetane, LHV, PMI, cost)</li>
                <li>• Enforces fuel-type specific constraints</li>
                <li>• Respects component caps and novel limits</li>
                <li>• Produces identical results for the same seed</li>
              </ul>
            </div>
            <div>
              <h3 className="mb-2 font-semibold text-primary">Constraints</h3>
              <p className="text-muted-foreground">
                The system enforces realistic fuel constraints:
              </p>
              <ul className="mt-2 space-y-1 text-muted-foreground">
                <li>• Volatility limits (T10/T90, RVP)</li>
                <li>• Minimum octane/cetane requirements</li>
                <li>• Density bounds (diesel/jet)</li>
                <li>• Component maximum volume fractions</li>
                <li>• Novel component caps</li>
              </ul>
            </div>
            <div>
              <h3 className="mb-2 font-semibold text-primary">Output Format</h3>
              <p className="text-muted-foreground">
                Results include volume fractions, mL and grams per 1000mL recipes, complete property analysis, 
                cost breakdown, and full metadata in JSON format for reproducibility.
              </p>
            </div>
          </CardContent>
        </Card>

        {/* Workflow */}
        <Card>
          <CardHeader>
            <CardTitle className="text-2xl">Unified Workflow</CardTitle>
          </CardHeader>
          <CardContent className="space-y-4">
            <p className="text-muted-foreground">
              The system bridges molecular discovery to practical formulation:
            </p>
            <div className="space-y-3">
              <div className="flex items-start gap-3">
                <div className="flex h-8 w-8 shrink-0 items-center justify-center rounded-full bg-primary/10 text-primary font-semibold">
                  1
                </div>
                <div>
                  <h3 className="font-semibold">Generate</h3>
                  <p className="text-sm text-muted-foreground">
                    Molecular generator produces candidate SMILES strings from fragment libraries
                  </p>
                </div>
              </div>
              <div className="flex items-start gap-3">
                <div className="flex h-8 w-8 shrink-0 items-center justify-center rounded-full bg-primary/10 text-primary font-semibold">
                  2
                </div>
                <div>
                  <h3 className="font-semibold">Convert</h3>
                  <p className="text-sm text-muted-foreground">
                    Conversion module estimates fuel properties and creates component records
                  </p>
                </div>
              </div>
              <div className="flex items-start gap-3">
                <div className="flex h-8 w-8 shrink-0 items-center justify-center rounded-full bg-primary/10 text-primary font-semibold">
                  3
                </div>
                <div>
                  <h3 className="font-semibold">Optimize</h3>
                  <p className="text-sm text-muted-foreground">
                    Blend optimizer selects and proportions components into spec-compliant blends
                  </p>
                </div>
              </div>
              <div className="flex items-start gap-3">
                <div className="flex h-8 w-8 shrink-0 items-center justify-center rounded-full bg-primary/10 text-primary font-semibold">
                  4
                </div>
                <div>
                  <h3 className="font-semibold">Validate</h3>
                  <p className="text-sm text-muted-foreground">
                    Outputs are immediately actionable for laboratory mixing and testing
                  </p>
                </div>
              </div>
            </div>
          </CardContent>
        </Card>
      </div>
    </div>
  );
}


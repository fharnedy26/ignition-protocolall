import Link from "next/link";
import { ArrowRight, Play, Code } from "lucide-react";
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card";

export default function DemoPage() {
  return (
    <div className="container mx-auto px-4 py-12">
      <div className="mb-12">
        <h1 className="mb-4 text-5xl font-bold">
          <span className="bg-gradient-to-r from-primary to-accent bg-clip-text text-transparent">
            Interactive Demo
          </span>
        </h1>
        <p className="text-xl text-muted-foreground">
          Try the fuel blend optimization engine with an interactive playground.
        </p>
      </div>

      <div className="mb-12 grid gap-6 md:grid-cols-2">
        <Link href="/demo/playground">
          <Card className="group h-full transition-all hover:border-primary hover:shadow-[0_0_20px_oklch(0.65_0.20_35/0.3)]">
            <CardHeader>
              <div className="mb-2 flex h-12 w-12 items-center justify-center rounded-lg bg-primary/10 text-primary">
                <Play className="h-6 w-6" />
              </div>
              <CardTitle>Playground</CardTitle>
              <CardDescription>
                Interactive interface to optimize fuel blends with component input, fuel type selection, and real-time results.
              </CardDescription>
            </CardHeader>
            <CardContent>
              <div className="flex items-center gap-2 text-primary group-hover:gap-4 transition-all">
                <span className="text-sm font-medium">Launch playground</span>
                <ArrowRight className="h-4 w-4" />
              </div>
            </CardContent>
          </Card>
        </Link>
        <Card>
          <CardHeader>
            <div className="mb-2 flex h-12 w-12 items-center justify-center rounded-lg bg-primary/10 text-primary">
              <Code className="h-6 w-6" />
            </div>
            <CardTitle>API Integration</CardTitle>
            <CardDescription>
              The playground uses the engine API to perform real optimizations.
            </CardDescription>
          </CardHeader>
          <CardContent>
            <p className="text-sm text-muted-foreground">
              The demo connects to the optimization engine via API routes, allowing you to input component data, 
              select fuel types, and receive optimized blend results with visualizations.
            </p>
          </CardContent>
        </Card>
      </div>

      <section>
        <h2 className="mb-4 text-3xl font-bold">Features</h2>
        <div className="grid gap-4 md:grid-cols-3">
          <Card>
            <CardHeader>
              <CardTitle className="text-lg">Component Input</CardTitle>
            </CardHeader>
            <CardContent>
              <p className="text-sm text-muted-foreground">
                Upload CSV files or enter component data manually with properties like density, RON, MON, cetane, and cost.
              </p>
            </CardContent>
          </Card>
          <Card>
            <CardHeader>
              <CardTitle className="text-lg">Fuel Type Selection</CardTitle>
            </CardHeader>
            <CardContent>
              <p className="text-sm text-muted-foreground">
                Choose from gasoline, diesel, or jet fuel specifications. Each type enforces appropriate constraints.
              </p>
            </CardContent>
          </Card>
          <Card>
            <CardHeader>
              <CardTitle className="text-lg">Visualizations</CardTitle>
            </CardHeader>
            <CardContent>
              <p className="text-sm text-muted-foreground">
                View blend composition, property comparisons, and component fractions with interactive charts.
              </p>
            </CardContent>
          </Card>
        </div>
      </section>
    </div>
  );
}


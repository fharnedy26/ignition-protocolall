"use client";

import { useState } from "react";
import { Button } from "@/components/ui/button";
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card";
import { Label } from "@/components/ui/label";
import { Input } from "@/components/ui/input";
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from "@/components/ui/select";
import { Textarea } from "@/components/ui/textarea";
import { Loader2, Upload, Play } from "lucide-react";
import { BlendComposition, PropertyComparison } from "@/components/charts";
import { motion, AnimatePresence } from "framer-motion";

export default function PlaygroundPage() {
  const [fuelType, setFuelType] = useState("gasoline");
  const [k, setK] = useState("5");
  const [top, setTop] = useState("10");
  const [components, setComponents] = useState("");
  const [loading, setLoading] = useState(false);
  const [results, setResults] = useState<any>(null);
  const [error, setError] = useState<string | null>(null);

  const handleOptimize = async () => {
    setLoading(true);
    setError(null);
    setResults(null);

    try {
      // Parse components CSV
      const lines = components.trim().split("\n");
      if (lines.length < 2) {
        throw new Error("Please provide at least one component");
      }

      // Simple CSV parsing (header + data)
      const headers = lines[0].split(",").map((h) => h.trim());
      const componentData = lines.slice(1).map((line) => {
        const values = line.split(",").map((v) => v.trim());
        const obj: any = {};
        headers.forEach((header, idx) => {
          obj[header] = values[idx] || "";
        });
        return obj;
      });

      const response = await fetch("/api/optimize", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({
          components: componentData,
          fuel_type: fuelType,
          K: parseInt(k),
          top: parseInt(top),
        }),
      });

      if (!response.ok) {
        const errorData = await response.json();
        throw new Error(errorData.error || "Optimization failed");
      }

      const data = await response.json();
      setResults(data);
    } catch (err: any) {
      setError(err.message || "An error occurred");
    } finally {
      setLoading(false);
    }
  };

  const exampleCSV = `id,name,density_g_ml,LHV_MJ_kg,RON,MON,cetane,PMI,TSI,cost_eur_L,max_vol_frac,is_novel
A,Gasoline,0.75,44.0,87,82,0,15,15,1.0,1.0,0
B,Ethanol,0.79,26.8,108,89,0,5,5,1.2,0.15,0
C,MTBE,0.74,35.2,118,101,0,10,10,1.5,0.20,0`;

  return (
    <div className="container mx-auto px-4 py-12">
      <motion.div
        initial={{ opacity: 0, y: 20 }}
        animate={{ opacity: 1, y: 0 }}
        transition={{ duration: 0.6 }}
        className="mb-12"
      >
        <h1 className="mb-4 text-5xl font-bold">
          <span className="bg-gradient-to-r from-primary to-accent bg-clip-text text-transparent">
            Optimization Playground
          </span>
        </h1>
        <p className="text-xl text-muted-foreground">
          Input component data and optimize fuel blends in real-time.
        </p>
      </motion.div>

      <div className="grid gap-8 lg:grid-cols-2">
        {/* Input Panel */}
        <motion.div
          initial={{ opacity: 0, x: -30 }}
          animate={{ opacity: 1, x: 0 }}
          transition={{ duration: 0.6, delay: 0.2 }}
          className="space-y-6"
        >
          <Card>
            <CardHeader>
              <CardTitle>Parameters</CardTitle>
              <CardDescription>Configure optimization settings</CardDescription>
            </CardHeader>
            <CardContent className="space-y-4">
              <div>
                <Label htmlFor="fuel-type">Fuel Type</Label>
                <Select value={fuelType} onValueChange={setFuelType}>
                  <SelectTrigger id="fuel-type">
                    <SelectValue />
                  </SelectTrigger>
                  <SelectContent>
                    <SelectItem value="gasoline">Gasoline</SelectItem>
                    <SelectItem value="diesel">Diesel</SelectItem>
                    <SelectItem value="jet">Jet Fuel</SelectItem>
                  </SelectContent>
                </Select>
              </div>
              <div>
                <Label htmlFor="k">Max Components (K)</Label>
                <Input
                  id="k"
                  type="number"
                  value={k}
                  onChange={(e) => setK(e.target.value)}
                  min="1"
                  max="10"
                />
              </div>
              <div>
                <Label htmlFor="top">Top Results</Label>
                <Input
                  id="top"
                  type="number"
                  value={top}
                  onChange={(e) => setTop(e.target.value)}
                  min="1"
                  max="20"
                />
              </div>
            </CardContent>
          </Card>

          <Card>
            <CardHeader>
              <CardTitle>Component Data</CardTitle>
              <CardDescription>Paste CSV data or use the example below</CardDescription>
            </CardHeader>
            <CardContent className="space-y-4">
              <div>
                <Label htmlFor="components">CSV Components</Label>
                <Textarea
                  id="components"
                  value={components}
                  onChange={(e) => setComponents(e.target.value)}
                  placeholder={exampleCSV}
                  className="min-h-[300px] font-mono text-sm"
                />
              </div>
              <Button
                variant="outline"
                onClick={() => setComponents(exampleCSV)}
                className="w-full"
              >
                <Upload className="mr-2 h-4 w-4" />
                Load Example
              </Button>
              <Button
                onClick={handleOptimize}
                disabled={loading || !components.trim()}
                className="w-full"
              >
                {loading ? (
                  <>
                    <Loader2 className="mr-2 h-4 w-4 animate-spin" />
                    Optimizing...
                  </>
                ) : (
                  <>
                    <Play className="mr-2 h-4 w-4" />
                    Optimize
                  </>
                )}
              </Button>
            </CardContent>
          </Card>
        </motion.div>

        {/* Results Panel */}
        <motion.div
          initial={{ opacity: 0, x: 30 }}
          animate={{ opacity: 1, x: 0 }}
          transition={{ duration: 0.6, delay: 0.4 }}
          className="space-y-6"
        >
          <AnimatePresence mode="wait">
            {error && (
              <motion.div
                key="error"
                initial={{ opacity: 0, scale: 0.9 }}
                animate={{ opacity: 1, scale: 1 }}
                exit={{ opacity: 0, scale: 0.9 }}
                transition={{ duration: 0.3 }}
              >
                <Card className="border-destructive">
                  <CardHeader>
                    <CardTitle className="text-destructive">Error</CardTitle>
                  </CardHeader>
                  <CardContent>
                    <p className="text-sm text-muted-foreground">{error}</p>
                  </CardContent>
                </Card>
              </motion.div>
            )}

            {results && (
              <motion.div
                key="results"
                initial={{ opacity: 0, y: 20 }}
                animate={{ opacity: 1, y: 0 }}
                exit={{ opacity: 0, y: -20 }}
                transition={{ duration: 0.5 }}
              >
            <Card>
              <CardHeader>
                <CardTitle>Optimization Results</CardTitle>
                <CardDescription>Best blend found</CardDescription>
              </CardHeader>
              <CardContent className="space-y-4">
                <div>
                  <Label className="text-primary">Score</Label>
                  <p className="text-2xl font-bold">{results.score?.toFixed(4) || "N/A"}</p>
                </div>
                <div>
                  <Label className="text-primary">Properties</Label>
                  <div className="mt-2 space-y-1 text-sm">
                    {results.properties && Object.entries(results.properties).map(([key, value]: [string, any]) => (
                      <div key={key} className="flex justify-between">
                        <span className="text-muted-foreground">{key}:</span>
                        <span className="font-mono">{typeof value === "number" ? value.toFixed(2) : value}</span>
                      </div>
                    ))}
                  </div>
                </div>
                <div>
                  <Label className="text-primary">Components</Label>
                  <div className="mt-2 space-y-1 text-sm">
                    {results.components?.map((comp: any, idx: number) => (
                      <div key={idx} className="flex justify-between">
                        <span className="text-muted-foreground">{comp.id || comp.name}:</span>
                        <span className="font-mono text-primary">
                          {(comp.vol_frac * 100).toFixed(1)}%
                        </span>
                      </div>
                    ))}
                  </div>
                </div>
              </CardContent>
            </Card>
              </motion.div>
            )}

            {results && results.components && (
              <motion.div
                key="visualizations"
                initial={{ opacity: 0, y: 20 }}
                animate={{ opacity: 1, y: 0 }}
                transition={{ duration: 0.5, delay: 0.2 }}
              >
                <Card>
                  <CardHeader>
                    <CardTitle>Visualizations</CardTitle>
                  </CardHeader>
                  <CardContent className="space-y-6">
                    <motion.div
                      initial={{ opacity: 0, scale: 0.9 }}
                      animate={{ opacity: 1, scale: 1 }}
                      transition={{ duration: 0.5, delay: 0.4 }}
                    >
                      <Label className="text-primary mb-2 block">Blend Composition</Label>
                      <BlendComposition components={results.components} />
                    </motion.div>
                    {results.properties && (
                      <motion.div
                        initial={{ opacity: 0, scale: 0.9 }}
                        animate={{ opacity: 1, scale: 1 }}
                        transition={{ duration: 0.5, delay: 0.6 }}
                      >
                        <Label className="text-primary mb-2 block">Property Comparison</Label>
                        <PropertyComparison properties={results.properties} />
                      </motion.div>
                    )}
                  </CardContent>
                </Card>
              </motion.div>
            )}

            {!results && !error && !loading && (
              <motion.div
                key="empty"
                initial={{ opacity: 0 }}
                animate={{ opacity: 1 }}
                exit={{ opacity: 0 }}
                transition={{ duration: 0.3 }}
              >
                <Card>
                  <CardContent className="pt-6">
                    <p className="text-center text-muted-foreground">
                      Configure parameters and input component data, then click Optimize to see results.
                    </p>
                  </CardContent>
                </Card>
              </motion.div>
            )}
          </AnimatePresence>
        </motion.div>
      </div>
    </div>
  );
}

